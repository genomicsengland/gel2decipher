import logging
from protocols.cva_1_0_0 import *

from protocols.migration.migration import Migration2_1To3
from protocols.migration.migration_reports_3_0_0_to_reports_4_0_0 import MigrateReports3To4
from protocols.migration.migration_reports_4_0_0_to_reports_5_0_0 import MigrateReports400To500

# Models 2.1.0 imports
from protocols.reports_2_1_0 import ClinicalReportRD as ClinicalReportRD_2_1_0
from protocols.reports_2_1_0 import InterpretedGenomeRD as InterpretedGenomeRD_2_1_0
from protocols.reports_2_1_0 import ClinicalReportCancer as ClinicalReportCancer_2_1_0
from protocols.reports_2_1_0 import InterpretationRequestRD as InterpretationRequestRD_2_1_0
# Models 3.0.0 imports
from protocols.reports_3_0_0 import ClinicalReportRD as ClinicalReportRD_3_0_0
from protocols.reports_3_0_0 import InterpretedGenomeRD as InterpretedGenomeRD_3_0_0
from protocols.reports_3_0_0 import ClinicalReportCancer as ClinicalReportCancer_3_0_0
from protocols.reports_3_0_0 import InterpretationRequestRD as InterpretationRequestRD_3_0_0
from protocols.reports_3_0_0 import CancerInterpretationRequest as CancerInterpretationRequest_3_0_0
# Models 4.0.0 imports
from protocols.reports_4_0_0 import ClinicalReportCancer as ClinicalReportCancer_4_0_0
from protocols.reports_4_0_0 import ClinicalReportRD as ClinicalReportRD_4_0_0
from protocols.reports_4_0_0 import InterpretationRequestRD as InterpretationRequestRD_4_0_0
from protocols.reports_4_0_0 import InterpretedGenomeRD as InterpretedGenomeRD_4_0_0
from protocols.reports_4_0_0 import CancerInterpretedGenome as CancerInterpretedGenome_4_0_0
from protocols.reports_4_0_0 import CancerInterpretationRequest as CancerInterpretationRequest_4_0_0
# Models 5.0.0
from protocols.reports_5_0_0 import Assembly as Assembly_5_0_0
from protocols.reports_5_0_0 import ClinicalReportRD as ClinicalReportRD_5_0_0
from protocols.reports_5_0_0 import InterpretedGenomeRD as InterpretedGenomeRD_5_0_0
from protocols.reports_5_0_0 import ClinicalReportCancer as ClinicalReportCancer_5_0_0
from protocols.reports_5_0_0 import InterpretationRequestRD as InterpretationRequestRD_5_0_0
from protocols.reports_5_0_0 import CancerInterpretedGenome as CancerInterpretedGenome_5_0_0
from protocols.reports_5_0_0 import CancerInterpretationRequest as CancerInterpretationRequest_5_0_0

from gel2decipher.clients.rest_client import RestClient
from gel2decipher.clients.model_validator import PayloadValidation


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
    def get_interpreted_genome(ir_result):
        ig_list = ir_result.get('interpreted_genome')
        sorted_ig_list = sorted(ig_list, key=lambda ig: ig['created_at'], reverse=True)
        latest_ig = next((ig for ig in sorted_ig_list), None)
        return latest_ig

    @staticmethod
    def get_clinical_report(ir_result):
        cr_list = ir_result.get('clinical_report')
        sorted_cr_list = sorted(cr_list, key=lambda cr: cr['created_at'], reverse=True)
        latest_cr = next((cr for cr in sorted_cr_list if cr != []), None)
        return latest_cr

    @staticmethod
    def split_assembly_from_patch(interpretation_request_rd_or_cancer):
        return next((assembly for assembly in interpretation_request_rd_or_cancer.genomeAssemblyVersion.split('.')),
                    None)

    def generate_tiered_variant_inject_rd(self, interpretation_request):
        """
        Tiered variants injection fields are described here:
        https://cnfl.extge.co.uk/display/IPI/CVA+mapping+data+from+the+CIPAPI+model+to+the+CVA+model
        """
        json_dict = interpretation_request['interpretation_request_data']['json_request']
        assembly = interpretation_request['assembly']
        interpreted_genome_rd = self.migrate_interpretation_request_rd_to_version_5_0_0(
            json_dict=json_dict, assembly=assembly)

        tiered_variant_inject_rd = TieredVariantInjectRD(
            id=str(interpretation_request.get('interpretation_request_id', None)),
            version=int(interpretation_request.get('version', None)),

            assembly=assembly,
            interpretedGenome=interpreted_genome_rd,
            author='tiering',
            # TODO: where does this field come from?
            authorVersion=getattr(interpreted_genome_rd, "tieringVersion", None),
            # TODO: where does this field come from?
            workspace=json_dict.get("workspace", None),  ## workspace is lost in the new InterpretedGenome
            groupId=interpretation_request.get('family_id', None),
            cohortId=interpretation_request.get('cohort_id', None),
            reportModelVersion=interpreted_genome_rd.versionControl.gitVersionControl,
        )
        return tiered_variant_inject_rd

    def generate_candidate_variant_inject_rd(self, interpreted_genome, interpretation_request):
        ir_json_dict = interpretation_request['interpretation_request_data']['json_request']
        assembly = interpretation_request['assembly']
        interpretation_request_rd = self.migrate_interpretation_request_rd_to_version_5_0_0(
            json_dict=ir_json_dict, assembly=assembly)
        version_control = interpretation_request_rd.versionControl.gitVersionControl
        workspace = ir_json_dict.get("workspace", None)  ## workspace is lost in the new InterpretedGenome
        parent_id = interpretation_request.get('interpretation_request_id', None)
        parent_version = interpretation_request.get('version', None)
        cohort_id = interpretation_request.get('cohort_id', None)
        group_id = interpretation_request.get('family_id', None)
        author = interpretation_request.get('cip', None)
        return self.generate_candidate_variant_inject_rd(
            interpreted_genome, assembly, version_control, workspace, parent_id, parent_version, cohort_id, group_id,
            author
        )

    def generate_candidate_variant_inject_rd(self, interpreted_genome, assembly, version_control, workspace,
                                             parent_id, parent_version, cohort_id, group_id, author):
        """
        Reported variants injection fields are described here:
        https://cnfl.extge.co.uk/display/IPI/CVA+mapping+data+from+the+CIPAPI+model+to+the+CVA+model
        """
        ig_json_dict = interpreted_genome['interpreted_genome_data']
        interpreted_genome_rd = self.migrate_interpreted_genome_rd_to_version_5_0_0(
            json_dict=ig_json_dict, assembly=assembly, interpretation_request_version=parent_version)

        ig_version = interpreted_genome['cip_version']
        candidate_variant_inject_rd = CandidateVariantInjectRD(
            # The id must "ir_id-ir_version"
            id="{parent_id}-{parent_version}".format(parent_id=parent_id, parent_version=parent_version),
            version=int(ig_version),

            assembly=assembly,
            interpretedGenome=interpreted_genome_rd,
            author=author,
            reportModelVersion=version_control,
            authorVersion=str(interpreted_genome['cip_version']),
            workspace=workspace,
            parentId=str(parent_id),
            cohortId=cohort_id,
            parentVersion=int(parent_version),
            groupId=group_id
        )
        return candidate_variant_inject_rd

    def generate_reported_variant_inject_rd(self, clinical_report, interpretation_request):
        """
        Candidate variants injection fields are described here:
        https://cnfl.extge.co.uk/display/IPI/CVA+mapping+data+from+the+CIPAPI+model+to+the+CVA+model
        """
        json_dict = interpretation_request['interpretation_request_data']['json_request']
        assembly = interpretation_request['assembly']
        interpretation_request_rd = self.migrate_interpretation_request_rd_to_version_5_0_0(
            json_dict=json_dict, assembly=assembly)

        version_control = interpretation_request_rd.versionControl.gitVersionControl
        workspace = json_dict.get("workspace", None)  ## workspace is lost in the new InterpretedGenome
        parent_id = interpretation_request.get('interpretation_request_id', None)
        parent_version = interpretation_request.get('version', None)
        cohort_id = interpretation_request.get('cohort_id', None)
        group_id = interpretation_request.get('family_id', None)

        return self.generate_reported_variant_inject_rd(
            clinical_report, assembly, version_control, workspace, parent_id, parent_version, cohort_id, group_id)

    def generate_reported_variant_inject_rd(self, clinical_report, assembly, version_control, workspace,
                                            parent_id, parent_version, cohort_id, group_id):
        """
        Candidate variants injection fields are described here:
        https://cnfl.extge.co.uk/display/IPI/CVA+mapping+data+from+the+CIPAPI+model+to+the+CVA+model
        """
        json_dict = clinical_report['clinical_report_data']
        clinical_report_rd = self.migrate_clinical_report_rd_to_version_5_0_0(json_dict=json_dict, assembly=assembly)

        cr_version = clinical_report['clinical_report_version']
        reported_variant_inject_rd = ReportedVariantInjectRD(
            # The id must "ir_id-ir_version"
            id="{parent_id}-{parent_version}".format(parent_id=parent_id, parent_version=parent_version),
            version=cr_version,
            clinicalReport=clinical_report_rd,
            assembly=assembly,
            author=getattr(clinical_report_rd, "user", None),
            parentId=str(parent_id),
            parentVersion=int(parent_version),
            reportModelVersion=version_control,
            authorVersion=None,
            workspace=workspace,
            groupId=group_id,
            cohortId=cohort_id,
        )
        return reported_variant_inject_rd

    @staticmethod
    def migrate_interpretation_request_rd_to_version_5_0_0(json_dict, assembly):
        ir_v500 = None

        if PayloadValidation(klass=InterpretedGenomeRD_5_0_0, payload=json_dict).is_valid:
            ir_v500 = InterpretedGenomeRD_5_0_0.fromJsonDict(jsonDict=json_dict)
            logging.info("Case in models reports 5.0.0")

        elif PayloadValidation(klass=InterpretationRequestRD_4_0_0, payload=json_dict).is_valid:
            ir_v400 = InterpretationRequestRD_4_0_0.fromJsonDict(jsonDict=json_dict)
            ir_v500 = MigrateReports400To500().migrate_interpretation_request_rd_to_interpreted_genome_rd(
                old_instance=ir_v400, assembly=assembly, interpretation_service="tiering",
                reference_database_versions={}, software_versions={}
            )
            logging.info("Case in models reports 4.0.0")

        elif PayloadValidation(klass=InterpretationRequestRD_3_0_0, payload=json_dict).is_valid:
            ir_v3 = InterpretationRequestRD_3_0_0.fromJsonDict(jsonDict=json_dict)
            ir_v400 = MigrateReports3To4().migrate_interpretation_request_rd(old_instance=ir_v3)
            ir_v500 = MigrateReports400To500().migrate_interpretation_request_rd_to_interpreted_genome_rd(
                old_instance=ir_v400, assembly=assembly, interpretation_service="tiering",
                reference_database_versions={}, software_versions={}
            )
            logging.info("Case in models reports 3.0.0")

        elif PayloadValidation(klass=InterpretationRequestRD_2_1_0, payload=json_dict).is_valid:
            ir_v2 = InterpretationRequestRD_2_1_0.fromJsonDict(jsonDict=json_dict)
            ir_v3 = Migration2_1To3().migrate_interpretation_request(interpretation_request=ir_v2)
            ir_v400 = MigrateReports3To4().migrate_interpretation_request_rd(old_instance=ir_v3)
            ir_v500 = MigrateReports400To500().migrate_interpretation_request_rd_to_interpreted_genome_rd(
                old_instance=ir_v400, assembly=assembly, interpretation_service="tiering",
                reference_database_versions={}, software_versions={}
            )
            logging.info("Case in models reports 2.1.0")

        if ir_v500 is not None:
            return ir_v500

        raise ValueError("INTERPRETATION REQUEST RD IS NOT IN VERSIONS: [2.1.0, 3.0.0, 4.0.0, 5.0.0]")

    @staticmethod
    def migrate_interpreted_genome_rd_to_version_5_0_0(json_dict, assembly, interpretation_request_version):
        ig_v500 = None

        if PayloadValidation(klass=InterpretedGenomeRD_5_0_0, payload=json_dict).is_valid:
            ig_v500 = InterpretedGenomeRD_5_0_0.fromJsonDict(jsonDict=json_dict)
            logging.info("Case in models reports 5.0.0")

        if PayloadValidation(klass=InterpretedGenomeRD_4_0_0, payload=json_dict).is_valid:
            ig_v400 = InterpretedGenomeRD_4_0_0.fromJsonDict(jsonDict=json_dict)
            ig_v500 = MigrateReports400To500().migrate_interpreted_genome_rd(
                old_instance=ig_v400, assembly=assembly, interpretation_request_version=interpretation_request_version
            )
            logging.info("Case in models reports 4.0.0")

        elif PayloadValidation(klass=InterpretedGenomeRD_3_0_0, payload=json_dict).is_valid:
            ig_v3 = InterpretedGenomeRD_3_0_0.fromJsonDict(jsonDict=json_dict)
            ig_v400 = MigrateReports3To4().migrate_interpreted_genome_rd(old_instance=ig_v3)
            ig_v500 = MigrateReports400To500().migrate_interpreted_genome_rd(
                old_instance=ig_v400, assembly=assembly, interpretation_request_version=interpretation_request_version
            )
            logging.info("Case in models reports 3.0.0")

        elif PayloadValidation(klass=InterpretedGenomeRD_2_1_0, payload=json_dict).is_valid:
            ig_v2 = InterpretedGenomeRD_2_1_0.fromJsonDict(jsonDict=json_dict)
            ig_v3 = Migration2_1To3().migrate_interpreted_genome(interpreted_genome=ig_v2)
            ig_v400 = MigrateReports3To4().migrate_interpreted_genome_rd(old_instance=ig_v3)
            ig_v500 = MigrateReports400To500().migrate_interpreted_genome_rd(
                old_instance=ig_v400, assembly=assembly, interpretation_request_version=interpretation_request_version
            )
            logging.info("Case in models reports 2.1.0")

        if ig_v500 is not None:
            return ig_v500

        raise ValueError("INTERPRETED GENOME RD IS NOT IN VERSIONS: [2.1.0, 3.0.0, 4.2.0]")

    @staticmethod
    def migrate_clinical_report_rd_to_version_5_0_0(json_dict, assembly):
        cr_v500 = None

        if PayloadValidation(klass=ClinicalReportRD_5_0_0, payload=json_dict).is_valid:
            cr_v500 = ClinicalReportRD_5_0_0.fromJsonDict(jsonDict=json_dict)
            logging.info("Case in models reports 5.0.0")

        elif PayloadValidation(klass=ClinicalReportRD_4_0_0, payload=json_dict).is_valid:
            cr_v4 = ClinicalReportRD_4_0_0.fromJsonDict(jsonDict=json_dict)
            cr_v500 = MigrateReports400To500().migrate_clinical_report_rd(
                old_instance=cr_v4, assembly=assembly
            )
            logging.info("Case in models reports 4.0.0")

        elif PayloadValidation(klass=ClinicalReportRD_3_0_0, payload=json_dict).is_valid:
            cr_v3 = ClinicalReportRD_3_0_0.fromJsonDict(jsonDict=json_dict)
            cr_v4 = MigrateReports3To4().migrate_clinical_report_rd(
                old_clinical_report_rd=cr_v3
            )
            cr_v500 = MigrateReports400To500().migrate_clinical_report_rd(
                old_instance=cr_v4, assembly=assembly
            )
            logging.info("Case in models reports 3.0.0")

        elif PayloadValidation(klass=ClinicalReportRD_2_1_0, payload=json_dict).is_valid:
            cr_v2 = ClinicalReportRD_2_1_0.fromJsonDict(jsonDict=json_dict)
            cr_v3 = Migration2_1To3().migrate_clinical_report(clinical_report=cr_v2)
            cr_v4 = MigrateReports3To4().migrate_clinical_report_rd(
                old_clinical_report_rd=cr_v3
            )
            cr_v500 = MigrateReports400To500().migrate_clinical_report_rd(
                old_instance=cr_v4, assembly=assembly
            )
            logging.info("Case in models reports 2.1.0")

        if cr_v500 is not None:
            return cr_v500

        raise ValueError("CLINICAL REPORT RD IS NOT IN VERSIONS: [2.1.0, 3.0.0, 4.0.0, 5.0.0]")
