import logging

from pbcommand.models.report import Report, Attribute


def run_hello_world(conditions, output_report):
    a = Attribute("hello", value='world')
    report = Report('pbinternal_hello_world', attributes=[a])
    report.write_json(output_report)
    return 0
