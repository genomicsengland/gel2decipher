#!/env/python
import argparse
import logging


def main():
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(description='Loader from the CIPAPI to CVA')
    parser.add_argument('--cipapi-url', help='The URL for the CIPAPI')
    parser.add_argument('--cipapi-user', help='The user for the CIPAPI')
    parser.add_argument('--cipapi-password', help='The password for the CIPAPI')
    parser.add_argument('--cva-url', help='The URL for CVA')
    parser.add_argument('--cva-user', help='The user for the CVA')
    parser.add_argument('--cva-password', help='The password for the CVA')
    parser.add_argument('--authentication-url', help='The URL for the authentication service to connect to CVA')
    args = parser.parse_args()

    config = {
        "cipapi_url": args.cipapi_url,
        "cipapi_user": args.cipapi_user,
        "cipapi_password": args.cipapi_password,
        "cva_url": args.cva_url,
        "cva_user": args.cva_user,
        "cva_password": args.cva_password,
        "authentication_url": args.authentication_url
    }
    #loader = CIPAPI2CVA(config)
    #loader.load()


if __name__ == '__main__':
    main()
