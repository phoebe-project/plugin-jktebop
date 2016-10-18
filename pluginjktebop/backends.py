from phoebe.parameters import ParameterSet
from phoebe import u, c
from phoebe.backend.backends import _extract_from_bundle_by_dataset, _extract_from_bundle_by_time

import numpy as np


def jktebop(b, compute, times=[], **kwargs):
    """
    Use John Southworth's code (jktebop) to compute radial velocities
    and light curves.  The code is available here:

    http://www.astro.keele.ac.uk/jkt/codes/jktebop.html

    jktebop must be installed and available on the system in order to
    use this plugin.

    Please see the link above for a list of publications to cite when using this
    code.

    According to jktebop's website:

        JKTEBOP models the two components as biaxial spheroids for the
        calculation of the reflection and ellipsoidal effects,
        and as spheres for the eclipse shapes.

    Note that the wrapper around jktebop only uses its forward model.
    Jktebop also includes its own fitting kinds, including bootstrapping.
    Those capabilities cannot be accessed from PHOEBE.

    Parameters that are used by this backend:

    * Compute:
        -  all parameters in :func:`phoebe.parameters.compute.jktebop`

    * Orbit:
        - TODO: list these

    * Star:
        - TODO: list these

    * lc dataset:
        - TODO :list these

    Values that are filled by this backend:

    * lc:
        - times
        - fluxes

    * rv (dynamical only):
        - times
        - rvs

    This function will almost always be called through the bundle, using
        * :meth:`phoebe.frontend.bundle.Bundle.add_compute`
        * :meth:`phoebe.frontend.bundle.Bundle.run_compute`

    :parameter b: the :class:`phoebe.frontend.bundle.Bundle` containing the system
        and datasets
    :parameter str compute: the label of the computeoptions to use (in the bundle).
        These computeoptions must have a kind of 'jktebop'.
    :parameter **kwargs: any temporary overrides to computeoptions
    :return: a list of new synthetic :class:`phoebe.parameters.parameters.ParameterSet`s
    :raise ImportError: if the jktebop executable cannot be found or is not installed
    :raises ValueError: if an ld_func is not valid for the jktebop backedn
    """

    # check whether jktebop is installed
    out = commands.getoutput('jktebop')
    if 'not found' in out:
        raise ImportError('jktebop executable not found')

    computeparams = b.get_compute(compute, force_ps=True)
    hier = b.get_hierarchy()

    starrefs  = hier.get_stars()
    orbitrefs = hier.get_orbits()

    if len(starrefs) != 2 or len(orbitrefs) != 1:
        raise ValueError("jktebop backend only accepts binary systems")

    logger.warning("JKTEBOP backend is still in development/testing and is VERY experimental")

    orbitref = orbitrefs[0]

    infos, new_syns = _extract_from_bundle_by_dataset(b, compute=compute, times=times)  # TODO: add rv support (see commented context below)

    ringsize = computeparams.get_value('ringsize', unit=u.deg, **kwargs)

    rA = b.get_value('rpole', component=starrefs[0], context='component', unit=u.solRad)
    rB = b.get_value('rpole', component=starrefs[1], context='component', unit=u.solRad)
    sma = b.get_value('sma', component=orbitref, context='component', unit=u.solRad)
    incl = b.get_value('incl', component=orbitref, context='component', unit=u.deg)
    q = b.get_value('q', component=orbitref, context='component')
    ecosw = b.get_value('ecosw', component=orbitref, context='component')
    esinw = b.get_value('esinw', component=orbitref, context='component')

    gravbA = b.get_value('gravb_bol', component=starrefs[0], context='component')
    gravbB = b.get_value('gravb_bol', component=starrefs[1], context='component')


    period = b.get_value('period', component=orbitref, context='component', unit=u.d)
    t0_supconj = b.get_value('t0_supconj', component=orbitref, context='component', unit=u.d)

    for info in infos:
        info = info[0]

        # get dataset-dependent things that we need
        l3 = b.get_value('l3', dataset=info['dataset'], context='dataset')
        # TODO: need to sum up pblums of each component - so need to write a function which will use the phoebe2 backend
        # to compute pblums that are coupled (or they need to be computed as constraints - I guess we'll see how fast that function runs)
        try:
            pblum = sum([b.get_value('pblum', dataset=info['dataset'], component=starref, context='dataset') for starref in starrefs])  # TODO: supposed to be in mags?
        except:
            raise ValueError("jktebop backend currently only supports decoupled pblums (b.set_value_all('pblum_ref', 'self'))")

        logger.warning("pblum in jktebop is sum of pblums (per-component): {}".format(pblum))
        pblum = -2.5 * np.log10(pblum) + 0.0

        ldfuncA = b.get_value('ld_func', component=starrefs[0], dataset=info['dataset'], context='dataset')
        ldfuncB = b.get_value('ld_func', component=starrefs[1], dataset=info['dataset'], context='dataset')

        ldcoeffsA = b.get_value('ld_coeffs', component=starrefs[0], dataset=info['dataset'], context='dataset')
        ldcoeffsB = b.get_value('ld_coeffs', component=starrefs[1], dataset=info['dataset'], context='dataset')

        if len(ldcoeffsA) != 2:
            logger.warning("ld_coeffs not compatible with jktebop - setting to (0.5,0.5)")
            ldcoeffsA = (0.5,0.5)
        if len(ldcoeffsB) != 2:
            logger.warning("ld_coeffs not compatible with jktebop - setting to (0.5,0.5)")
            ldcoeffsB = (0.5,0.5)

        albA = b.get_value('frac_refl_bol', component=starrefs[0], context='component')
        albB = b.get_value('frac_refl_bol', component=starrefs[1], context='component')

        tratio = b.get_value('teff', component=starrefs[0], context='component', unit=u.K) / b.get_value('teff', component=starrefs[1], context='component', unit=u.K)

        # provide translation from phoebe's 'ld_func' to jktebop's 'LD law type'
        ldfuncs = {'linear': 'lin',
                'logarithmic': 'log',
                'square_root': 'sqrt',
                'quadratic': 'quad'}

        # let's make sure we'll be able to make the translation later
        if ldfuncA not in ldfuncs.keys() or ldfuncB not in ldfuncs.keys():
            # TODO: provide a more useful error statement
            raise ValueError("jktebop only accepts the following options for ld_func: {}".format(ldfuncs.keys()))

        # create the input file for jktebop
        fi = open('_tmp_jktebop_in', 'w')
        #~ fi.write("# JKTEBOP input file created by PHOEBE\n")

        # We always want task 2 - according to jktebop's website:
        # Task 2    This inputs a parameter file and calculates a
        # synthetic light curve (10000 points between phases 0 and 1)
        # using the parameters you put in the file.
        fi.write('{:5} {:11} Task to do (from 1 to 9)   Integ. ring size (deg)\n'.format(2, ringsize))
        fi.write('{:5} {:11} Sum of the radii           Ratio of the radii\n'.format((rA+rB)/sma, rA/rB))
        fi.write('{:5} {:11} Orbital inclination (deg)  Mass ratio of the system\n'.format(incl, q))

        # we'll provide ecosw and esinw instead of ecc and long_an
        # jktebop's readme.txt states that so long as ecc is < 10,
        # it will be intrepreted as ecosw and esinw (otherwise would
        # need to be ecc+10 and long_an (deg)
        fi.write('{:5} {:11} Orbital eccentricity       Periastron longitude deg\n'.format(ecosw, esinw))


        fi.write('{:5} {:11} Gravity darkening (starA)  Grav darkening (starB)\n'.format(gravbA, gravbB))
        fi.write('{:5} {:11} Surface brightness ratio   Amount of third light\n'.format(tratio, l3))


        # According to jktebop's readme.txt:
        # The possible entries for the type of limb darkening law are 'lin' (for linear)
        # 'log' (logarithmic), 'sqrt' (square-root), 'quad' (quadratic) or 'cub' (cubic)

        fi.write('{:5} {:11} LD law type for star A     LD law type for star B\n'.format(ldfuncs[ldfuncA], ldfuncs[ldfuncB]))
        fi.write('{:5} {:11} LD star A (linear coeff)   LD star B (linear coeff)\n'.format(ldcoeffsA[0], ldcoeffsB[0]))
        fi.write('{:5} {:11} LD star A (nonlin coeff)   LD star B (nonlin coeff)\n'.format(ldcoeffsA[1], ldcoeffsB[1]))

        fi.write('{:5} {:11} Reflection effect star A   Reflection effect star B\n'.format(albA, albB))
        fi.write('{:5} {:11} Phase of primary eclipse   Light scale factor (mag)\n'.format(0.0, pblum))
        fi.write('{:13}      Orbital period of eclipsing binary system (days)\n'.format(period))
        fi.write('{:13}      Reference time of primary minimum (HJD)\n'.format(t0_supconj))

        # All fitting will be done with PHOEBE wrappers, so we need to set
        # all jktebop options for adjust to False (0)
        fi.write(' {:d}  {:d}             Adjust RADII SUM or RADII RATIO (0, 1, 2, 3)\n'.format(0, 0))
        fi.write(' {:d}  {:d}             Adjust INCLINATION or MASSRATIO (0, 1, 2, 3)\n'.format(0, 0))
        fi.write(' {:d}  {:d}             Adjust ECCENTRICITY or OMEGA (0, 1, 2, 3)\n'.format(0, 0))
        fi.write(' {:d}  {:d}             Adjust GRAVDARK1 or GRAVDARK2 (0, 1, 2, 3)\n'.format(0, 0))
        fi.write(' {:d}  {:d}             Adjust SURFACEBRIGH2 or THIRDLIGHT (0, 1, 2, 3)\n'.format(0, 0))
        fi.write(' {:d}  {:d}             Adjust LD-lin1 or LD-lin2 (0, 1, 2, 3)\n'.format(0, 0))
        fi.write(' {:d}  {:d}             Adjust LD-nonlin1 or LD-nonlin2 (0, 1, 2, 3)\n'.format(0, 0))
        fi.write(' {:d}  {:d}             Adjust REFLECTION COEFFS 1 and 2 (-1, 0, 1, 2, 3)\n'.format(0, 0))
        fi.write(' {:d}  {:d}             Adjust PHASESHIFT or SCALE FACTOR (0, 1, 2, 3)\n'.format(0, 0))
        fi.write(' {:d}  {:d}             Adjust PERIOD or TZERO (min light) (0, 1, 2, 3)\n'.format(0, 0))

        #~ fi.write('{}  Name of file containing light curve\n'.format("_tmp_jktebop_lc_in"))
        #~ fi.write('{}  Name of output parameter file\n'.format("_tmp_jktebop_param"))
        #~ fi.write('{}  Name of output light curve file\n'.format("_tmp_jktebop_lc_out"))
        #~ fi.write('{}  Name of output model light curve fit file\n'.format("_tmp_jktebop_modelfit_out"))
        #~ fi.write('{}\n')

        # According to jktebop's readme.txt:
        # FITTING FOR RADIAL VELOCITIES:    the observed RVs should be in separate files
        # for the two stars and the data should be in the same format as the light curve
        # data. Then add a line below the main input parameters for each rv file:
        #   RV1  [infile]  [outfile]  [K]  [Vsys]  [vary(K)]  [vary(Vsys)]
        #   RV2  [infile]  [outfile]  [K]  [Vsys]  [vary(K)]  [vary(Vsys)]
        # where RV1 is for primary star velocities, RV2 is for secondary star velocities
        # [infile] is the input data file, [outfile] is the output data file, [K] is the
        # velocity amplitude of the star (km/s), [Vsys] is its systemic velocity (km/s),
        # and [vary(K)] and [vary(Vsys)] are 0 to fix and 1 to fit for these quantities.
        # The mass ratio parameter is not used for the RVs, only for the light curve.
        # If you want to fix the systemic velocity for star B to that for star A, simply
        # set vary(Vsys) for star B to be equal to -1
        #~ fi.write('rv1 llaqr-rv1.dat llaqr-rv1.out 55.0 -10.0 1 1\n')
        #~ fi.write('rv2 llaqr-rv2.dat llaqr-rv2.out 55.0 -10.0 1 1\n')


        # According to jktebop's readme.txt:
        # NUMERICAL INTEGRATION:  long exposure times can be split up into NUMINT points
        # occupying a total time interval of NINTERVAL (seconds) by including this line:
        #   NUMI  [numint]  [ninterval]

        # TODO: allow exposure times


        fi.close()

        # TODO: create_tmp_jktebop_lc_in - probably with times and dummy fluxes if none are in the obs
        #~ flc = open('_tmp_jktebop_lc_in', 'w')
        #~ times = b.get_value('times', component=info['component'], dataset=info['dataset'], context='dataset', unit=u.d)
        #~ fluxes = b.get_value('flux', component=info['component'], dataset=info['dataset'], context='dataset', unit=u.d)

        #~ if len(fluxes) < len(times):
            #~ # then just provide dummy fluxes - we're not using
            #~ # jktebop for fitting anyways, it really just needs the times
            #~ fluxes = [1.]*len(times)

        #~ for t,f in zip(times, fluxes):
            #~ flc.write('{}\t{}\n'.format(t,t))
        #~ flc.close()

        # run jktebop
        out = commands.getoutput("jktebop _tmp_jktebop_in > _tmp_jktebop_lc_out")

        # parse output
        phases_all, mags_all, l1, l2, l3 = np.loadtxt(str(period), unpack=True)
        #~ time, flux = np.loadtxt("_tmp_jktebop_lc_out", unpack=True)

        # fill syn
        this_syn = new_syns.filter(component=info['component'], dataset=info['dataset'])

        # phases_all, mags_all are 10001 evenly-spaced phases, so we need to interpolate
        # to get at the desired times
        times_all = b.to_time(phases_all)  # in days
        mags_interp = np.interp(info['times'], times_all, mags_all)

        this_syn['times'] = info['times'] * u.d # (period was requested in days)
        logger.warning("converting from mags from JKTEBOP to flux")
        ref_mag = 0  # TODO: what should we do with this?? - option in jktebop compute?
        this_syn['fluxes'] = 10**((mags_interp-ref_mag)/-2.5) * 2  # 2 seems to be necessary - probably from a difference in pblum conventions (or just normalization???)

    yield new_syns
