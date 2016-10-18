from phoebe.parameters import *
from phoebe import u

def jktebop(**kwargs):
    """
    Compute options for using John Southworth's 'jktebop' code as a
    backend (must be installed).

    Generally, this will be used as an input to the kind argument in
    :meth:`phoebe.frontend.bundle.Bundle.add_compute`

    Please see :func:`phoebe.backend.backends.jktebop` for a list of sources to
    cite when using this backend.

    :parameter **kwargs: defaults for the values of any of the parameters
    :return: a :class:`phoebe.parameters.parameters.ParameterSet` of all newly
        created :class:`phoebe.parameters.parameters.Parameter`s
    """
    if not conf.devel:
        raise NotImplementedError("'jktebop' backend not officially supported for this release.  Enable developer mode to test.")

    params = []

    params += [BoolParameter(qualifier='enabled', copy_for={'context': 'dataset', 'kind': ['LC'], 'dataset': '*'}, visible_if='False', dataset='_default', value=kwargs.get('enabled', True), description='Whether to create synthetics in compute/fitting run')]

    params += [FloatParameter(qualifier='ringsize', value=kwargs.get('ringsize', 5), default_unit=u.deg, description='Integ Ring Size')]

    return ParameterSet(params)
