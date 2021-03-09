class attrdict(dict):
    """
    Dictionary which provides attribute access to its keys.
    """

    def __getattr__(self, key):
        if key in self:
            return self[key]
        else:
            raise AttributeError(
                "%r object has no attribute %r" % (type(self).__name__, key)
            )
