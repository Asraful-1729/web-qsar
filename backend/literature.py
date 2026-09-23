"""
A4 — Literature intelligence: fetches PubMed papers (title/abstract/authors)
relevant to a free-text query (a target name, compound name, or plant
source) via NCBI's E-utilities (esearch + efetch), per explicit user
choice of PubMed over other literature APIs.

This is the app's FIRST live, request-time external network dependency —
A2's curated pool is bundled and A3's COCONUT index is downloaded once
ahead of time, but this hits the network on every search. Every call is
therefore defensive about connectivity failures / timeouts and raises a
distinct LiteratureError with a plain-language message (mapped to HTTP
502 in app.py) rather than a raw exception, since a desktop app user may
simply be offline.

No API key is used: NCBI's documented rate limit without one is 3
requests/second, comfortably above what a single-user desktop app
searching a handful of times per session needs. See
https://www.ncbi.nlm.nih.gov/books/NBK25497/ for the E-utilities contract.
"""
import xml.etree.ElementTree as ET

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
TIMEOUT = 12  # seconds -- fail fast rather than hang a request-time UI call


class LiteratureError(Exception):
    """Network/API failure talking to PubMed — distinct from ValueError
       (bad input) so app.py can map it to a 502, not a 400."""


def _get(client, url, params):
    import httpx  # local import matches this project's existing lazy-import convention for optional deps
    try:
        r = client.get(url, params=params, timeout=TIMEOUT)
        r.raise_for_status()
        return r.content
    except httpx.TimeoutException as e:
        raise LiteratureError("PubMed request timed out.") from e
    except httpx.HTTPStatusError as e:
        raise LiteratureError(f"PubMed returned an error ({e.response.status_code}).") from e
    except httpx.HTTPError as e:
        raise LiteratureError(f"Could not reach PubMed ({e}). Check your internet connection.") from e


def search(query, max_results=10):
    """Returns {query, n_results, papers}, or raises ValueError (empty
       query) / LiteratureError (network/API failure)."""
    import httpx
    import json as _json

    if not query or not query.strip():
        raise ValueError("empty query")
    query = query.strip()
    max_results = max(1, min(int(max_results), 30))

    with httpx.Client() as client:
        raw = _get(client, ESEARCH_URL, {
            "db": "pubmed", "term": query, "retmax": max_results, "retmode": "json", "sort": "relevance",
        })
        try:
            pmids = _json.loads(raw)["esearchresult"]["idlist"]
        except (KeyError, ValueError) as e:
            raise LiteratureError("Unexpected response from PubMed search.") from e

        if not pmids:
            return {"query": query, "n_results": 0, "papers": []}

        xml_raw = _get(client, EFETCH_URL, {
            "db": "pubmed", "id": ",".join(pmids), "rettype": "abstract", "retmode": "xml",
        })

    try:
        root = ET.fromstring(xml_raw)
    except ET.ParseError as e:
        raise LiteratureError("Could not parse PubMed response.") from e

    papers = []
    for article in root.findall(".//PubmedArticle"):
        pmid_el = article.find(".//PMID")
        pmid = pmid_el.text if pmid_el is not None else None

        title_el = article.find(".//ArticleTitle")
        title = "".join(title_el.itertext()).strip() if title_el is not None else "(no title)"

        abstract_parts = article.findall(".//AbstractText")
        abstract = " ".join("".join(p.itertext()).strip() for p in abstract_parts) if abstract_parts else None

        journal_el = article.find(".//Journal/Title")
        journal = journal_el.text if journal_el is not None else None

        year_el = article.find(".//PubDate/Year")
        if year_el is None:
            year_el = article.find(".//PubDate/MedlineDate")
        year = year_el.text[:4] if year_el is not None and year_el.text else None

        all_authors = article.findall(".//AuthorList/Author")
        authors = []
        for a in all_authors[:3]:
            last = a.find("LastName")
            if last is not None and last.text:
                authors.append(last.text)
        author_str = ", ".join(authors)
        if len(all_authors) > 3:
            author_str += " et al."

        papers.append({
            "pmid": pmid,
            "title": title,
            "abstract": abstract,
            "journal": journal,
            "year": year,
            "authors": author_str,
            "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else None,
        })

    # efetch's XML order isn't guaranteed to match esearch's relevance
    # ranking, so re-sort by esearch's original idlist order.
    order = {pmid: i for i, pmid in enumerate(pmids)}
    papers.sort(key=lambda p: order.get(p["pmid"], len(pmids)))

    return {"query": query, "n_results": len(papers), "papers": papers}
