#returns current date/time string in MMM-DD-YYYY HH:MM:SS.SSSSSS
import datetime
def vms_date():
	return  datetime.datetime.now().strftime("%b-%d-%Y %H:%M:%S.%f").upper()[:-4]
