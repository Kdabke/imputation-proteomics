#!/bin/bash
qsub -N Docker /common/dabkek/openswath_newlib/docker/docker.job
qsub -hold_jid Docker -N OpenSwath /common/dabkek/openswath_newlib/openswath/openswath.job
qsub -hold_jid OpenSwath -N PyProphet /common/dabkek/openswath_newlib/pyprophet/pyprophet.job
qsub -hold_jid PyProphet -N FeatureAlignment /common/dabkek/openswath_newlib/featurealignment/featurealignment.job
