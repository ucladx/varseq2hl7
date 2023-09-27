from subprocess import run

# given a list of file paths to VarSeq JSONs, send them to the flask server that sends HL7 messages to beaker

with open("hl7_jsons.txt", "r") as f:
	jsons = f.readlines()
	for j in jsons:
		j = j.strip()
		print(f"Sending {j}")
		run(["curl", "-s", "-o", "/dev/null", "-d", "@"+j, "-H", "Content-Type: application/json", "-X", "POST", "http://localhost:5000/receivejson"])
		cont = input("Enter to continue, anything else to stop: ")
		if cont != "":
			break
