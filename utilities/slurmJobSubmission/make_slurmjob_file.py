import sys
import argparse


def make_new_line(line, args):
    new_line = line
    for key, value in args.__dict__.items():
        if value is not None:
            new_line = new_line.replace("{"+key+"}", value)

    return new_line


def main(args):

    slurmjob_template_directory = "slurmjob_templates"
    slurmjob_template_file = slurmjob_template_directory + "/" + args.slurmJobTemplateName + ".sh"

    with open(slurmjob_template_file, "r") as file_:
        template = file_.readlines()

    with open(args.fileName, "w") as slurmjob_file:
        for line in template:
            new_line = make_new_line(line, args)
            slurmjob_file.write(new_line)

    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s", "--slurmJobTemplateName",
        )
    parser.add_argument(
        "-p", "--partition",
        )
    parser.add_argument(
        "-t", "--timeLimit",
        )
    parser.add_argument(
        "-m", "--memory",
        )
    parser.add_argument(
        "-l", "--logDir",
        )
    parser.add_argument(
        "-i", "--inputFile",
        )
    parser.add_argument(
        "-o", "--outputFile",
        )
    parser.add_argument(
        "-f", "--fileName",
        )
    parser.add_argument(
        "-j", "--jobName",
        )
    parser.add_argument(
        "-xrdcp", "--xrdcpInputFile",
        default="false"
        )
    parser.add_argument(
        "-sample", "--sample",
        )
    parser.add_argument(
        "-xsec", "--xsec",
        )


    main(parser.parse_args())
