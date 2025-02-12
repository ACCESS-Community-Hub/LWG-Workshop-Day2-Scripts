import os


def list_nc_files(d):
    """Lists all .nc files in subdirectories, or as files themselves."""
    if os.path.isfile(d):
        return [d]

    files = []
    for r, d, f in os.walk(d):
        for file in f:
            if ".nc" in file:
                files.append(os.path.join(r, file))
    return files
