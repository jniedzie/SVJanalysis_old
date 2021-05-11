import sys

def make_new_line(line, argDict):
    newLine = line
    for key in argDict.keys():
        newLine = newLine.replace("{"+key+"}", argDict[key])

    return newLine


def main(sample, argDict, outFileName):

    slurmjobTemplateFile = "slurmjob_template_" + sample + ".sh"

    with open(slurmjobTemplateFile, "r") as f:
        template = f.readlines()

    with open(outFileName, "w") as f:
        for line in template:
            newLine = make_new_line(line, argDict)
            f.write(newLine)

    return 0


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Incorrect usage!")
        sys.exit(1)

    if sys.argv[1] == "QCD":

        sample = "QCD"

        if len(sys.argv) != 6:
            print("Incorrect usage!")
            sys.exit(1)

        argDict = {
            "Pt1": sys.argv[2],
            "Pt2": sys.argv[3],
            "dataset": sys.argv[4],
        }

        outFileName = sys.argv[5]


    elif sys.argv[1] == "QCD_schannelCuts":

        sample = "QCD_schannelCuts"

        if len(sys.argv) != 7:
            print("Incorrect usage!")
            sys.exit(1)

        argDict = {
            "pt1": sys.argv[2],
            "pt2": sys.argv[3],
            "f1": sys.argv[4],
            "f2": sys.argv[5],
        }

        outFileName = sys.argv[6]


    if sys.argv[1] == "tchannel" or sys.argv[1] == "schannel":

        sample = sys.argv[1]

        if len(sys.argv) != 7:
            print("Incorrect usage!")
            sys.exit(1)

        argDict = {
            "mMed": sys.argv[2],
            "mDark": sys.argv[3],
            "rinv": sys.argv[4],
            "alpha": sys.argv[5],
        }

        outFileName = sys.argv[6]


    main(sample, argDict, outFileName)
