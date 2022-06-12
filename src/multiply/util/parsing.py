import os
import configparser
import datetime
from multiply.download.collection import genome_collection
from multiply.util.exceptions import DesignFileError


def parse_parameters(design_path):
    """
    Parse command-line input parameters passed to MULTIPLY

    params
        design_path: str
            Path to MULTIPLY design file.


    """

    config = configparser.ConfigParser()
    config.read(design_path)

    check_design_exists(design_path)
    check_valid_sections(config)

    params = {}
    params = add_samples(config, params)
    params = add_genes(config, params)
    params = add_regions(config, params)
    params = add_primers(config, params)
    params = add_amplicons(config, params)
    params = add_outputs(config, params)

    check_genes_or_regions(params)

    return params


def check_design_exists(design_path):
    """
    Check that the `design_path` points to an `.ini` file

    params
        design_path: str
            Path to MULTIPLY design file.
    returns
        None

    """

    if not os.path.isfile(design_path):
        raise FileNotFoundError(
            f"No design file found at {design_path}. Presently at {os.getcwd()}. Check the path."
        )
    if not design_path.endswith(".ini"):
        raise DesignFileError(f"The design {design_path} should end with `.ini`.")


def check_valid_sections(
    config,
    must_include=["Sample", "Primers", "Amplicons", "Output"],
    one_of=[["Genes", "Regions"]],
):
    """
    Check that all expected sections are found within the configuration object

    params
        config: ConfigParser
            ConfigParser object holding design file information.
        must_include: list
            List of sections that must be included in the `config` object.
        one_of: list of lists
            Each sublist gives a set of sections, at least one of which
            must be included in `config` object.

    returns
        None

    """

    for section in must_include:
        if not config.has_section(section):
            raise DesignFileError(f"Design missing the [{section}] section. Please add.")

    for section_set in one_of:
        has_one = False
        for section in section_set:
            if config.has_section(section):
                has_one = True
        if not has_one:
            raise DesignFileError(
                f"Design must include at least one of these sections: {', '.join(section_set)}. Please add."
            )


def add_samples(config, params):
    """
    Add [Sample] information to a parameter dictionary

    params
        config: ConfigParser
            ConfigParser object holding design file information.
        params: dict
            Dictionary of MULTIPLY parameters.
    returns
        params: dict
            Dictionary of MULTIPLY parameters, with [Sample]
            parameters added.

    """

    # Check that it is within the enumerated samples
    genome = config.get("Sample", "genome")

    if genome not in genome_collection:
        raise DesignFileError(f"In [Sample], the provided `genome` {genome} is not in collection. Run 'multiply download --available' to see available genomes.")
        # Must be in the registry

    params["genome"] = genome

    return params


def add_genes(config, params):
    """
    Add [Genes] information from a configparser object to a params dictionary

    params
        config: ConfigParser
            ConfigParser object holding design file information.
        params: dict
            Dictionary of MULTIPLY parameters.
    returns
        params: dict
            Dictionary of MULTIPLY parameters, with [Genes]
            parameters added.

    """

    # Check if genes have been provided
    if not config.has_section("Genes"):
        params["from_genes"] = False
        return params

    # Parse gene IDs
    target_ids = [g.strip() for g in config.get("Genes", "target_ids").split(",")]

    # Parse gene namess
    has_names = config.has_option("Genes", "target_names")
    if has_names:
        target_names = [
            g.strip() for g in config.get("Genes", "target_names").split(",")
        ]
    else:
        target_names = target_ids

    # Mapping between IDs and names
    id_to_name = {i: n for i, n in zip(target_ids, target_names)}
    name_to_id = {n: i for i, n in id_to_name.items()}

    # Sanity checks
    n_ids = len(target_ids)
    n_names = len(target_names)
    if not n_ids == n_names:
        raise DesignFileError(
            f"In [Genes], found {n_ids} `target_ids` and {n_names} `target_names`. Ensure equal."
        )

    # Add to dictionary
    params["from_genes"] = True
    params["target_ids"] = target_ids
    params["has_names"] = has_names
    params["target_names"] = target_names
    params["target_id_to_name"] = id_to_name

    return params


def add_regions(config, params):
    """
    Add [Regions] information from a configparser object to a params dictionary

    params
        config: ConfigParser
            ConfigParser object holding design file information.
        params: dict
            Dictionary of MULTIPLY parameters.
    returns
        params: dict
            Dictionary of MULTIPLY parameters, with [Regions]
            parameters added.

    """

    # Check if regions have been provided
    if not config.has_section("Regions"):
        params["from_regions"] = False
        return params

    params["from_regions"] = True
    params["region_bed"] = config.get("Regions", "bed")

    return params


def add_primers(config, params):
    """
    Add [Primer] information from a configparser object to a params dictionary

    """

    params["include_tails"] = config.getboolean("Primers", "include_tails")
    if params["include_tails"]:
        
        # TODO: probably add sanity check here? Or no?
        params["F_tail"] = config.get("Primers", "F_tail").upper()
        params["R_tail"] = config.get("Primers", "R_tail").upper()

    return params


def add_amplicons(config, params, min_size_bp=50, max_size_bp=10000):
    """
    Add [Amplicon] section to a parameter dictionary

    params
        config: ConfigParser
            ConfigParser object holding design file information.
        params: dict
            Dictionary of MULTIPLY parameters.
    returns
        params: dict
            Dictionary of MULTIPLY parameters, with [Amplicons]
            parameters added.

    """

    params["min_size_bp"] = config.getint("Amplicons", "min_size_bp")
    if params["min_size_bp"] < min_size_bp:
        print(
            f"In [Amplicon], min_size_bp = {params['min_size_bp']} is too small. Setting to {min_size_bp}bp."
        )
        params["min_size_bp"] = min_size_bp
    params["max_size_bp"] = config.getint("Amplicons", "max_size_bp")
    if params["max_size_bp"] > max_size_bp:
        print(
            f"In [Amplicon], max_size_bp = {params['max_size_bp']} is too large. Setting to {max_size_bp}bp."
        )
        params["max_size_bp"] = max_size_bp

    # primer3 settings
    if config.has_option("Amplicons", "primer3_settings"):
        params["primer3_settings"] = [
            s.strip() for s in config.get("Amplicons", "primer3_settings").split(",")
        ]
    else:
        params["primer3_settings"] = ["default"]

    return params


def add_outputs(config, params):
    """
    Add [Output] section to a parameter dictionary

    params
        config: ConfigParser
            ConfigParser object holding design file information.
        params: dict
            Dictionary of MULTIPLY parameters.
    returns
        params: dict
            Dictionary of MULTIPLY parameters, with [Output]
            parameters added.

    """

    params["output_name"] = config.get("Output", "name")
    if config.has_option("Output", "primer_code"):
        params["primer_code"] = config.get("Output", "primer_code")
    else:
        params["primer_code"] = "m"

    # Create an output directory
    today = datetime.datetime.today().strftime("%Y-%m-%d")
    params["output_dir"] = f"results/{today}_{params['output_name']}"

    return params


def check_genes_or_regions(params):
    """
    Check that a parameter set has either [Genes] or [Regions]
    defined

    """

    if not (params["from_genes"] or params["from_regions"]):
        raise DesignFileError(f"Either [Genes] and/or [Regions] must be specified.")
