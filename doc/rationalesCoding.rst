------------------------
Coding-style conventions
------------------------

#. Reserve ``struct`` for TMP (traits classes, metafunctions, etc.) and stateless classes, otherwise always use ``class``

#. Names representing types must be in mixed case starting with upper case

#. Variable names must be in mixed case starting with lower case

#. Private class variables must have underscore suffix (this makes them easier to spot than the underscore prefix.)

#. Named constants (including enumeration values, and template parameters of integral const type) must be all uppercase using underscore to separate words --- Enumeration constants must be prefixed by a common type name

#. Names representing methods or functions should be verbs and written in mixed case starting with lower case --- They are distinguished from variable names anyway by their special syntax

#. Names representing template metafunctions must be named as if they were functions, but in mixed case starting with upper case as if they were types

#. Names representing namespaces should be all lowercase

#. Names representing template types should be a single uppercase letter

#. Generic variables should have the same name as their type --- Non-generic variables have a role: these variables can often be named by combining role and type

#. Variables with a large scope should have long names, variables with a small scope can have short names

#. The name of the object is implicit, and should be avoided in a method name

#. The terms `get/set` must be used where an attribute is accessed directly --- This, together with items 2 3 4 6 10 allows for very convenient naming, eg::
  
	class Foo
	{
	public:
	  Foo(Bar bar) : bar_(bar) {}  

	  bar  getBar(       ) const {return bar_;}
	  void setBar(Bar bar)       {bar_=bar   ;}
	private:
	  Bar bar_;
	};




#. Plural form should be used on names representing a collection of objects

#. The prefix *n* should be used for variables representing a number of objects

#. The suffix *No* should be used for variables representing an entity number

#. Iterator variables should be called i, j, k etc.

#. The prefix *is* should be used for boolean variables and methods

#. Naming pointers specifically should be avoided

#. Negated boolean variable names must be avoided

#. Exception classes should be suffixed with ``Exception``

#. Functions (methods returning something) should be named after what they return, procedures (void methods) after what they do

#. C++ header files should have the extension ``.h``, source files can have the extension ``.cpp`` or ``.cc``

#. Header files must contain an include guard

#. Include statements should be sorted and grouped.  Sorted by their hierarchical position in the system with high level files included first. Leave an empty line between groups of include statements. Cf. below in code organization

#. Include statements must be located at the top of a file only

#. The parts of a class must be sorted ``public``, ``protected``, and ``private``.  All sections must be identified explicitly. Not applicable sections should be left out.

#. Abbreviations and acronyms must not be uppercase when used as name

#. Complement names must be used for complement operations

#. C++ pointers and references should have their reference symbol next to the type rather than to the name

#. The nominal case should be put in the if-part and the exception in the else-part of an if statement

#. Use alignment wherever it enhances readability

#. Use // for all comments, including multi-line comments --- then any larger passage can be commented out with ``/*...*/``

#. The function return type can be put in the left column immediately above the function name

#. Consider using a beautifier (bcpp?)

.. #. The conditional should be put on a separate line
