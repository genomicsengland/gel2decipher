#!/env/python
import argparse
import logging

from gel2decipher.case_sender import Gel2Decipher


def main():
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(description='Loader from the CIPAPI to CVA')
    parser.add_argument('--cipapi-url', help='The URL for the CIPAPI', required=True)
    parser.add_argument('--cva-url', help='The URL for CVA', required=True)
    parser.add_argument('--gel-user', help='The user for GEL', required=True)
    parser.add_argument('--gel-password', help='The password for GEL', required=True)
    parser.add_argument('--decipher-system-key', help="Decipher's system key", required=True)
    parser.add_argument('--decipher-user-key', help="Decipher's user key", required=True)
    parser.add_argument('--decipher-url', help="Decipher's URL", required=True)
    parser.add_argument('--send-absent-phenotypes', help="Flag to send absent phenotypes", action='store_true')
    args = parser.parse_args()

    config = {
        "cipapi_url": args.cipapi_url,
        "cva_url": args.cva_url,
        "gel_user": args.gel_user,
        "gel_password": args.gel_password,
        "decipher_system_key": args.decipher_system_key,
        "decipher_user_key": args.decipher_user_key,
        "decipher_url": args.decipher_url,
        "send_absent_phenotypes": args.send_absent_phenotypes
    }
    loader = Gel2Decipher(config)
    loader.send_case()


if __name__ == '__main__':
    main()
