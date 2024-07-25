import gspread
from google.oauth2.service_account import Credentials
scopes = ["https://www.googleapis.com/auth/spreadsheets"]
creds  = Credentials.from_service_account_file("kn1d-project-credentials-spreadsheets-135c1abb3187.json",scopes=scopes)
client = gspread.authorize(creds)

sheet_id = "1bdl3KibDaG6Yxah7OL7B39nl0ReJws950qINjtHt84g"
sheet = client.open_by_key(sheet_id)

values_list = sheet.sheet1.row_values(1)
print(values_list)