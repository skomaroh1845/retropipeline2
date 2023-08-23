import requests

API_link = "https://api.telegram.org/bot"

updates = requests.get(API_link + "/getUpdates").json()

print(updates)