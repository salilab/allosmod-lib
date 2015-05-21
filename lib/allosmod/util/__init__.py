def read_templates(template_file):
    """Read list of templates (one per line) from the given file"""
    with open(template_file) as fh:
        return [line.rstrip('\r\n') for line in fh]

