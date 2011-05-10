*********
ChangeLog
*********

2010/05/10 Switching to a more refined definition of time-dependence schemes in :class:`structure::Hamiltonian`

  As a consequence, :class:`structure::TimeIndependentHamiltonian` has been deprecated.


2010/02/15 Switching to `Bazaar <https://sourceforge.net/scm/?type=bzr&group_id=187775>`_ revision control

  History has not been migrated to Bazaar, the old history remaining accessible from `CVS <https://sourceforge.net/scm/?type=cvs&group_id=187775>`_. Cf. the `tutorial <http://cppqed.sourceforge.net/tutorial/installation.html#obtaining-c-qed>`_ for further details.


2010/11/29 Change in the quantumtrajectory bundle

  Classes in the quantumtrajectory bundle (``MCWF_Trajectory``, ``Master``, etc.) do not have the IS_NO template parameter any more. This is because the selection between orthogonal and non-orthogonal state vectors / density operators will now be done relying on run-time polymorphy instead of compile-time polymorphy as so far.


2010/08/26 Change in ParameterTable

  A change occured in how parameters::ParameterTable handles boolean parameters. For each parameter, it now automatically adds another parameter with a ``no_`` prefix which corresponds to the negation of the boolean. It is implemented in terms of a ``BooleanNegatedProxy`` which can be found in the header file of the same name. This makes that now there it makes no sense to declare negated boolean parameters any more. Therefore, e.g. the so far ubiquitous parameter ``nonoise`` (with default value ``false``) has been changed to ``noise`` (with default value ``true``) while the parameter ``no_noise`` is automatically added by ParameterTable (with default value ``false``, of course).


2010/08/12 Change in the API of Mode

  A change occured in the API of Mode as it got migrated to a template-based solution where the Averaging class is a plugin supplied as a template parameter. Furthermore, LossyModes acquired another template parameter signifying whether their temperature is finite. Both template parameters have default values, but still, at places where Modes are constructed explicitely (e.g. in scripts) the following change needs to be effected::

    PumpedLossyMode mode(...); // to be changed to
    PumpedLossyMode<> mode(...);

  where ``<>`` signifies that the class is a template, although with default template arguments. The use of the ``mode::maker`` function is unchanged::

    mode::SmartPtr mode(maker(...)); // is fine

  Similar changes will occur in Particle, Spin, Qbit, as they too get migrated to this solution.
