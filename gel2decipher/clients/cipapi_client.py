import logging

from gel2decipher.clients.rest_client import RestClient
from protocols.migration.migration_helpers import MigrationHelpers
from protocols.participant_1_0_3 import Pedigree, PedigreeMember


class CipApiClient(RestClient):

    ENDPOINT_BASE = "api/2"
    AUTH_ENDPOINT = "{}/get-token/".format(ENDPOINT_BASE)
    IR_ENDPOINT = "{}/interpretation-request/{}/{}"
    IR_LIST_ENDPOINT = "{}/interpretation-request?page={}&page_size=100"

    def __init__(self, url_base, token=None, user=None, password=None):
        """
        If user and password are not provided there will be no token renewal
        :param url_base:
        :param token:
        :param user:
        :param password:
        """
        RestClient.__init__(self, url_base)
        self.token = "JWT {}".format(token) if token is not None else None
        self.user = user
        self.password = password
        if (self.token is None) and (self.user is None or self.password is None):
            raise ValueError("Authentication is required. Provide either token or user and password.")
        self.set_authenticated_header()

    def get_token(self):
        token = self.post(self.AUTH_ENDPOINT, payload={
            'username': self.user,
            'password': self.password
        }).get('token')
        return "JWT {}".format(token)

    def get_interpretation_requests(self):
        page = 1
        while True:
            try:
                results = self.get(endpoint=self.IR_LIST_ENDPOINT.format(self.ENDPOINT_BASE, page))["results"]
            except ValueError:
                logging.info("Finished iterating through report events in the CIPAPI")
                break
            for result in results:
                last_status = result["last_status"]
                if last_status in ["blocked", "waiting_payload"]:
                    continue
                case_id, case_version = result["interpretation_request_id"].split("-")
                assembly = result["assembly"]
                yield self.get_interpretation_request(case_id, case_version)
            page += 1

    def get_interpretation_request(self, ir_id, ir_version):
        interpretation_request_result = self.get(
            endpoint=self.IR_ENDPOINT.format(self.ENDPOINT_BASE, ir_id, ir_version))

        return interpretation_request_result

    @staticmethod
    def get_interpreted_genome(ir_json):
        ig_list = ir_json.get('interpreted_genome')
        sorted_ig_list = sorted(ig_list, key=lambda ig: ig['created_at'], reverse=True)
        latest_ig = next((ig for ig in sorted_ig_list), None)
        return latest_ig

    @staticmethod
    def get_clinical_report(ir_json):
        cr_list = ir_json.get('clinical_report')
        sorted_cr_list = sorted(cr_list, key=lambda cr: cr['created_at'], reverse=True)
        latest_cr = next((cr for cr in sorted_cr_list if cr != []), None)
        return latest_cr

    @staticmethod
    def get_pedigree(ir_json):
        pedigree_json = ir_json['interpretation_request_data']['json_request']['pedigree']
        pedigree = MigrationHelpers.migrate_pedigree_to_version_1_0_3(pedigree_json)
        return pedigree

    @staticmethod
    def get_proband(pedigree):
        """

        :param pedigree:
        :type pedigree: Pedigree
        :rtype: PedigreeMember
        :return:
        """
        proband = None
        for participant in pedigree.participants:
            if participant.isProband:
                proband = participant
        return proband

    @staticmethod
    def split_assembly_from_patch(interpretation_request_rd_or_cancer):
        return next((assembly for assembly in interpretation_request_rd_or_cancer.genomeAssemblyVersion.split('.')),
                    None)
