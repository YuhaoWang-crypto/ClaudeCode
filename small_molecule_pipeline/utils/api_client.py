"""HTTP utilities: rate-limited requests with retry and exponential back-off."""

import time
import logging
from typing import Any, Dict, Optional
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

logger = logging.getLogger(__name__)

DEFAULT_TIMEOUT = 30
DEFAULT_PAUSE = 0.5  # seconds between requests (respect public APIs)


def _make_session(retries: int = 3, backoff: float = 1.0) -> requests.Session:
    session = requests.Session()
    retry = Retry(
        total=retries,
        backoff_factor=backoff,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET", "POST"],
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers.update({"User-Agent": "SmallMoleculePipeline/1.0 (research)"})
    return session


_SESSION = _make_session()


def get(
    url: str,
    params: Optional[Dict] = None,
    headers: Optional[Dict] = None,
    pause: float = DEFAULT_PAUSE,
    timeout: int = DEFAULT_TIMEOUT,
) -> Optional[requests.Response]:
    time.sleep(pause)
    try:
        resp = _SESSION.get(url, params=params, headers=headers, timeout=timeout)
        resp.raise_for_status()
        return resp
    except requests.HTTPError as e:
        logger.warning("HTTP error %s for %s: %s", e.response.status_code, url, e)
    except requests.RequestException as e:
        logger.warning("Request failed for %s: %s", url, e)
    return None


def post(
    url: str,
    data: Optional[Any] = None,
    json: Optional[Any] = None,
    headers: Optional[Dict] = None,
    pause: float = DEFAULT_PAUSE,
    timeout: int = DEFAULT_TIMEOUT,
) -> Optional[requests.Response]:
    time.sleep(pause)
    try:
        resp = _SESSION.post(url, data=data, json=json, headers=headers, timeout=timeout)
        resp.raise_for_status()
        return resp
    except requests.HTTPError as e:
        logger.warning("HTTP error %s for %s: %s", e.response.status_code, url, e)
    except requests.RequestException as e:
        logger.warning("Request failed for %s: %s", url, e)
    return None
