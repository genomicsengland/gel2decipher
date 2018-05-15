import logging
import requests
import datetime
import json
import abc
from requests.compat import urljoin
from requests.exceptions import HTTPError
import gel2decipher.clients.backoff_retrier as backoff_retrier


class RestClient(object):

    session = requests.Session()

    def __init__(self, url_base, retries=5):
        self.url_base = url_base
        self.headers = {
            'Accept': 'application/json'
        }
        self.token = None
        self.renewed_token = False
        # decorates the REST verbs with retries
        self.get = backoff_retrier.wrapper(self.get, retries)
        self.post = backoff_retrier.wrapper(self.post, retries)
        self.delete = backoff_retrier.wrapper(self.delete, retries)

    @staticmethod
    def build_url(baseurl, endpoint):
        return urljoin(baseurl, endpoint)

    def set_authenticated_header(self, renew_token=False):
        if not self.token or renew_token:
            self.token = self.get_token()
        self.headers["Authorization"] = "{token}".format(token=self.token)

    @abc.abstractmethod
    def get_token(self):
        raise ValueError("Not implemented")

    def post(self, endpoint, payload, url_params={}, session=True):
        if endpoint is None or payload is None:
            raise ValueError("Must define payload and endpoint before post")
        url = self.build_url(self.url_base, endpoint)
        logging.debug("{date} {method} {url}".format(
            date=datetime.datetime.now(),
            method="POST",
            url="{}?{}".format(url, "&".join(["{}={}".format(k, v) for k, v in url_params.iteritems()]))
        ))
        if session:
            response = self.session.post(url, json=payload, params=url_params, headers=self.headers)
        else:
            response = requests.post(url, json=payload, params=url_params, headers=self.headers)
        self._verify_response(response)
        return json.loads(response.content) if response.content else None

    def patch(self, endpoint, payload, url_params={}, session=True):
        if endpoint is None or payload is None:
            raise ValueError("Must define payload and endpoint before post")
        url = self.build_url(self.url_base, endpoint)
        logging.debug("{date} {method} {url}".format(
            date=datetime.datetime.now(),
            method="PATCH",
            url="{}?{}".format(url, "&".join(["{}={}".format(k, v) for k, v in url_params.iteritems()]))
        ))
        if session:
            response = self.session.patch(url, json=payload, params=url_params, headers=self.headers)
        else:
            response = requests.patch(url, json=payload, params=url_params, headers=self.headers)
        self._verify_response(response)
        return json.loads(response.content) if response.content else None

    def get(self, endpoint, url_params={}, session=True):
        if endpoint is None:
            raise ValueError("Must define endpoint before get")
        url = self.build_url(self.url_base, endpoint)
        logging.debug("{date} {method} {url}".format(
            date=datetime.datetime.now(),
            method="GET",
            url="{}?{}".format(url, "&".join(["{}={}".format(k, v) for k, v in url_params.iteritems()]))
        ))
        if session:
            response = self.session.get(url, params=url_params, headers=self.headers)
        else:
            response = requests.get(url, params=url_params, headers=self.headers)
        self._verify_response(response)
        return json.loads(response.content) if response.content else None

    def delete(self, endpoint, url_params={}):
        if endpoint is None:
            raise ValueError("Must define endpoint before get")
        url = self.build_url(self.url_base, endpoint)
        logging.debug("{date} {method} {url}".format(
            date=datetime.datetime.now(),
            method="DELETE",
            url="{}?{}".format(url, "&".join(["{}={}".format(k, v) for k, v in url_params.iteritems()]))
        ))
        response = self.session.delete(url, params=url_params, headers=self.headers)
        self._verify_response(response)
        return json.loads(response.content) if response.content else None

    def _verify_response(self, response):
        logging.debug("{date} response status code {status}".format(
            date=datetime.datetime.now(),
            status=response.status_code)
        )
        if response.status_code != 200:
            logging.error(response.content)
            # first 403 renews the token, second 403 in a row fails
            if response.status_code == 403 and not self.renewed_token:
                # renews the token if unauthorised
                self.set_authenticated_header(renew_token=True)
                self.renewed_token = True
                # RequestException will trigger a retry and with the renewed token it may work
                raise requests.exceptions.RequestException(response=response)
            # ValueError will not
            raise HTTPError("{}:{}".format(response.status_code, response.text), response=response)
        else:
            # once a 200 response token is not anymore just renewed, it can be renewed again if a 403 arrives
            self.renewed_token = False
