"""One-time CAVE auth: opens the browser to create a token, then saves it locally.
Accept the MICrONS terms first: https://global.daf-apis.com/sticky_auth/api/v1/tos/2/accept
Note: make_new=True invalidates any older token."""
from getpass import getpass
from caveclient import CAVEclient

client = CAVEclient()
client.auth.setup_token(make_new=True)
token = getpass("Paste the token from the browser (input hidden): ").strip()
client.auth.save_token(token=token, overwrite=True)
print("Token saved.")
