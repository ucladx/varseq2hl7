from subprocess import run
import argparse

# given a list of file paths to VarSeq JSONs (hl7_jsons.txt), send them to the flask server that sends HL7 messages to beaker
parser = argparse.ArgumentParser(description='Send VarSeq JSONs to Flask server')
parser.add_argument('--jsons_file', help='Path to the file containing the list of VarSeq JSONs')
args = parser.parse_args()

with open(args.jsons_file, "r") as f:
	jsons = f.readlines()
	for j in jsons:
		j = j.strip()
		print(f"Sending {j}")
		run(["curl", "-s", "-o", "/dev/null", "-d", "@"+j, "-H", "Content-Type: application/json", "-X", "POST", "http://localhost:5000/receivejson"])
		cont = input("Enter to continue, anything else to stop: ")
		if cont != "":
			break
