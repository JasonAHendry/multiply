class DesignFileError(Exception):
    """ Error in format or contents of the design file """
    pass

class NoTargetsFoundError(Exception):
    """ Error finding targets, usually related to GFF """
    pass

class TargetSizeError(Exception):
    """ Error in size of a target, usually relative to specified parameters """
    pass

class TargetPositionError(Exception):
    """ Error in the position of targets, i.e. some overlap """
    pass

class NoPrimerNameException(Exception):
    """ Raised when the `primer_name` column is missing """
    pass

class NoPrimersFoundException(Exception):
    """ Raised  when no primers are found for a target """
    pass
