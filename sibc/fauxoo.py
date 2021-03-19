"""
This shows the current state of my attempt at making a mixin class which allows
writing mixin-style code with a "global" (not really global) scope instead of
an explicit "self" instance.

It sets the global scope for each mixin suite using "exec" (of a code object).

Data and functions from multiple mixin suites can be combined. They share
read-write access to the same (instance-specific) global scope. It can also be
combined with normal OO classes.

>>> foobar = lambda: foo().use(bar)
>>> barfoo = lambda: bar().use(foo)
>>> foobar().collision()
'bar'
>>> barfoo().collision()
'foo'
>>> o1 = foobar()
>>> o2 = foobar()
>>> o1.add_AB()
3
>>> o1.B = 30
>>> o1.set_foriegn_A(7)
7
>>> o1.A
7
>>> o2.A
1
>>> o1.add_AB()
37
>>> o1.inc_A()
8
>>> o1.add_AB()
38

We can combine mixins with a normal class:

>>> n1 = normal()
>>> n1.set_foriegn_baz(23)
23
>>> n1.baz
23
>>> n1.get_baz()
23
>>> n1.add_AB()
3
>>> n1.inc_A()
2

Here we have a mixin of just one suite, and add the variable it wanted from
outside:

>>> a=bar()
>>> a.get_foriegn_A()
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "mixin1.py", line 138, in get_foriegn_A
    return A
NameError: name 'A' is not defined
>>> a.A=123
>>> a.get_foriegn_A()
123
"""

class fauxclass(object):

    def __init__(self, f):
        self._fn = f
        self.name = f.__name__

    def __call__(self):
        return fauxobj().use(self)

class fauxobj(object):

    """
    Faux objects are a construct similar to python's native classes, in that
    they are collections of methods and state which can be instantiated and
    composed. Unlike a normal Python class, however, a fauxclass is defined as
    function which is executed with an object instance as its global
    dictionary.

    This means that the 'self' variable is implicit: methods' function
    signatures do not include self, and attributes are accessed as global
    variables.

    When a method of a fauxclass needs to assign a value to a name in the
    instance dictionary, instead of using self.name, it can use the 'global'
    keyword, which will bind the name in the instance dictionary (and *not* the
    module global scope).

    At the end of the outer fauxclass definition function, after all methods
    have been defined, there should be a call to 'self.register(locals())'
    (this is the only use of self which is necessary within a fauxclass; self
    is however in fact always present with all of the correct properties.)
    """

    def __init__(self, _dict=None):
        if _dict is None:
            _dict = vars(self)
        self._fauxobj_dict = _dict
        self._fauxobj_bases = [ type(self).__name__ ]
        self.self = self # this makes self available to the fauxclasses as self

    def use(self, f):
        """
        This applies a fauxclass to a fauxobject.

        After this is called, the fauxobject could be said to inherit from the
        fauxclass.
        """
        assert isinstance(f, fauxclass)
        self._fauxobj_bases.append(f.name)
        # this is what makes it all work: we exec the fauxclass function with
        # our vars(self) as its globals() - this means that its globals
        # dictionary is actually the *very same dictionary* as our vars(self)
        # (aka self.__dict__)
        exec(f._fn.__code__, self._fauxobj_dict)
        return self

    def register(self, methods):
        """
        This merges the methods (typically locals(), but could be another
        dictionary with a subset of the methods) of the fauxclass into the
        instance dictionary.

        This should be called at the end of every fauxclass definition.
        """
        for k, v in methods.items():
            if not callable(v):
                raise ValueError('Refusing to register a non-callable: %s=%r' % (k,v))
        self._fauxobj_dict.update(methods)

    def __repr__(self):
        data = self._fauxobj_dict.copy()
        data.pop('__builtins__',None)
        data.pop('_fauxobj_bases',None)
        del data['self']
        methods = []
        for k,v in list(data.items()):
            if callable(v):
                methods.append(k)
                del data[k]
        return "<mix of %r with methods %s and attributes %s>" % (self._fauxobj_bases, methods, data)

    @classmethod
    def monkeypatch(cls, obj, f):
        """
        Bonus feature: This lets us apply fauxclasses to already-existing
        non-participating objects. One caveat is that when applied this way,
        the optional 'self' global in the fauxclass will not have the object's
        instance attributes on it (but will still be there to call register
        on).

        Use with caution.
        """
        vars(obj).setdefault('self', cls(vars(obj))).use(f)


if __name__ == "__main__":

    @fauxclass
    def foo():
        global A
        A = 1
        def add(v):
            return v + A
        def inc_A():
            global A
            A+=1
            return A
        def set_new(v):
            global AAA
            AAA = v
        def set_new_using_self(v):
            # this is equivalent, and also still works
            self.BBB = v
        def get_foriegn_B():
            return B
        def collision():
            return "foo"
        self.register(locals())

    @fauxclass
    def bar():
        global B
        B = 2
        def add_B(v):
            return v + B
        def inc_B():
            global B
            B+=1
            return B
        def add_AB():
            return A + B
        def get_foriegn_A():
            return A
        def set_foriegn_A(v):
            global A
            A = v
            return A
        def get_foriegn_baz():
            return baz
        def set_foriegn_baz(v):
            global baz
            baz = v
            return baz
        def collision():
            return "bar"
        self.register(locals())

    class normal(fauxobj):
        def __init__(self):
            self.baz = 23
            super(normal,self).__init__()
            self.use(foo).use(bar)

        def inc_both(self):
            self.inc_A()
            self.inc_B()
        
        def get_baz(self):
            return self.baz
        
        def set_baz(self, v):
            self.baz = v

        def get_A_normal(self):
            return self.A
        def collision(self):
            return "normal"

    class unwitting(object):
        def __init__(self):
            self.baz = 23

        def inc_both(self):
            self.inc_A()
            self.inc_B()
        
        def get_baz(self):
            return self.baz
        
        def set_baz(self, v):
            self.baz = v

        def get_A_normal(self):
            return self.A
        def collision(self):
            return "unwitting"

    import doctest
    doctest.testmod()
