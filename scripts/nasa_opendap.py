"""Helpers for spatially subsetting NASA Earthdata granules through Hyrax OPeNDAP (DAP4).

Authentication uses earthaccess with strategy="netrc" (credentials file created by the
user themselves — this code never asks for or stores passwords).
"""
from __future__ import annotations

import time
import xml.etree.ElementTree as ET
from urllib.parse import quote

import earthaccess

NS = {"d": "http://xml.opendap.org/ns/DAP/4.0#"}


def login():
    auth = earthaccess.login(strategy="netrc")
    if not auth or not auth.authenticated:
        raise SystemExit(
            "Earthdata login failed: no credentials file found.\n"
            "Run once in your own terminal (you will be prompted; nothing is stored in this project):\n"
            "  .venv\\Scripts\\python -c \"import earthaccess; earthaccess.login(strategy='interactive', persist=True)\""
        )
    return earthaccess.get_requests_https_session()


def opendap_url(granule) -> str:
    """Cloud Hyrax URL built from the collection concept id + full GranuleUR.

    The 'OPENDAP DATA' RelatedUrl stored in CMR for GPM_3IMERGHH V07 is truncated
    (drops '.0000.V07B.HDF5') and returns 404, so we do not use it.
    """
    from urllib.parse import quote as _q
    cid = granule["meta"]["collection-concept-id"]
    return f"https://opendap.earthdata.nasa.gov/collections/{cid}/granules/{_q(granule['umm']['GranuleUR'], safe='')}"


def get(session, url, retries=4, **kw):
    for i in range(retries):
        try:
            r = session.get(url, timeout=180, **kw)
            r.raise_for_status()
            return r
        except Exception:
            if i == retries - 1:
                raise
            time.sleep(2 ** i * 3)


def var_dims(session, base_url: str) -> dict[str, list[str]]:
    """Parse the DMR: {'/Group/var': ['/dim1', '/dim2', ...]}."""
    root = ET.fromstring(get(session, base_url + ".dmr.xml").content)
    out = {}

    def walk(node, prefix):
        for child in node:
            tag = child.tag.split("}")[-1]
            name = child.get("name")
            if tag == "Group":
                walk(child, f"{prefix}/{name}")
            elif name and child.findall("d:Dim", NS):
                out[f"{prefix}/{name}"] = [d.get("name") for d in child.findall("d:Dim", NS)]
            elif name and tag in ("Float32", "Float64", "Int32", "Int16", "Int64", "UInt8", "Int8", "UInt16"):
                out.setdefault(f"{prefix}/{name}", [])

    walk(root, "")
    return out


def subset_url(base_url: str, ce_parts: list[str]) -> str:
    return f"{base_url}.dap.nc4?dap4.ce=" + quote(";".join(ce_parts), safe="")
