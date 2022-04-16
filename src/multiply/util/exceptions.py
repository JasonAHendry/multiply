class DesignFileError(Exception):
    """ Error in format or contents of the design file """
    pass

class TargetSizeError(Exception):
    """ Error in size of a target, usually relative to specified parameters """
    pass

class TargetPositionError(Exception):
    """ Error in the position of targets, i.e. some overlap """
    pass
