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
    slurmjob_template_file = slurmjob_template_directory + "/" + args.sample + ".sh"

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
        "-s", "--sample",
        default="default",
        )
    parser.add_argument(
        "-l", "--logDir",
        )
    parser.add_argument(
        "-i", "--inputFile",
        )
    parser.add_argument(
        "-xsec", "--xsec",
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


    main(parser.parse_args())
