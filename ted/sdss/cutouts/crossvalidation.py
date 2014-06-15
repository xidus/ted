#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version: Fri 28 Mar 2014
#   Initial build.
#

import os
import time

import numpy as np
import pandas as pd
import yaml

# ... : ted
# from ... import msg
from ... import env
# from ... import TEDError
from ... import Namespace

# . : ted.sdss.cutouts
from . import load_cutout_sequences


class CVHelper(object):
    """Loads parameters and makes various parameters and methods available"""

    def __init__(self):
        self._load_parameters()
        self._load_data()
        self._set_cs_parameters()
        self.N_folds = self._cv.get('N_folds')

    # Initialise
    # ----------

    def _load_parameters(self):
        self._params = load_parameters()

    def _load_data(self):
        self._css, self._targets = load_cutout_sequences()

    def _set_cs_parameters(self):
        """Set attribute values for bg_model, clip,
        load_quality and quality for each CS."""
        for cs in self._css:
            cs.set_cs_parameters(**self._cs)

    @property
    def _cv(self): return self._params.get('cv')

    @property
    def _cs(self): return self._params.get('cs')

    @property
    def _xp(self): return self._params.get('xp')

    # Experiment properties
    # ---------------------

    def set_exp(self, key=None):
        if self._xp.get(key) is None: return None
    #     getattr(self, '_set_exp_' + key)()

    # def _set_exp_any(self):
    #     xpp = dict2numpy(self._xp.get('any'))
        # xpp contains the parameter grid alone.
        xpp = dict2numpy(self._xp.get(key))
        xp = Namespace()
        xp.name = key
        xp.sigmas = xpp.get('sigmas')
        xp.taus = xpp.get('taus')
        xp.N_sigmas = xp.sigmas.size
        xp.N_taus = xp.taus.size
        xp.gskw = dict(exp=xp.name, **xpp)
        self.xp = xp

    # QUALITY

    @property
    def quality(self):
        if not hasattr(self, '_quality'):
            # Use the one specified in the parameter file
            return self._cs.get('quality')
        else:
            # If it is set alreadt, use it.
            return self._quality

    def set_quality(self, quality):
        if isinstance(quality, list):
            self._quality = quality

    @staticmethod
    def _qstr(quality):
        return ''.join([str(q) for q in sorted(quality)])

    @property
    def qstr(self):
        return self._qstr(self.quality)

    # CUTOUT-SEQUENCE FILENAMES

    @property
    def _fname_prediction(self):
        # return 'prediction_Q-{}.csv'.format(self.qstr)
        return 'prediction-E-{}_Q-{}.csv'.format(self.xp.name, self.qstr)

    # @property
    def _fname_signals(self, sigma, tau):
        fstr = 'signals-E-{}_Q-{}_S-{:.2f}_T-{:.2f}.csv'
        return fstr.format(self.xp.name, self.qstr, sigma, tau)

    # Fold management
    # ---------------

    @property
    def _folds(self): return get_folds(self._css, self._targets, self.N_folds)

    # Cross validation
    # ----------------

    def cv(self):
        self._cv_time = time.time()
        getattr(self, '_cv_' + self.xp.name)()

    # ANY

    def _cv_any(self):
        for n, data in enumerate(self._folds):

            # Fold number
            fnum = n + 1

            # Unpack training and test fold (f: features; l: labels)
            (train_f, train_l), (test_f, test_l) = data

            try:
                # Training set
                cop = self._cop_any(train_f)
                moa = self._moa(cop, train_l)
                self._save_fold_results_moa(moa=moa, ftype='train', fnum=fnum)

            except Exception as e:
                print '_cv_any:', fnum, 'train'
                print '_cv_any:', e.message
                print '_cv_any:', '-' * 50, '\n'

            try:
                # Test set
                cop = self._cop_any(test_f)
                moa = self._moa(cop, test_l)
                self._save_fold_results_moa(moa=moa, ftype='test', fnum=fnum)

            except Exception as e:
                print '_cv_any:', fnum, 'test'
                print '_cv_any:', e.message
                print '_cv_any:', '-' * 50, '\n'

    # ANY2

    def _cv_any2(self): self._cv_any()

    # BASELINE

    def _cv_blr(self): self._cv_baseline()
    def _cv_bla(self): self._cv_baseline()
    def _cv_bln(self): self._cv_baseline()

    def _cv_blr2(self): self._cv_baseline()
    def _cv_bla2(self): self._cv_baseline()
    def _cv_bln2(self): self._cv_baseline()

    def _cv_baseline(self):

        for n, data in enumerate(self._folds):

            # Fold number
            fnum = n + 1

            # Unpack training and test fold (f: features; l: labels)
            (train_f, train_l), (test_f, test_l) = data

            # Training set
            cop = self._cop_baseline(train_f)
            moa = self._moa(cop, train_l)
            self._save_fold_results_moa(moa=moa, ftype='train', fnum=fnum)

            # Test set
            cop = self._cop_baseline(test_f)
            moa = self._moa(cop, test_l)
            self._save_fold_results_moa(moa=moa, ftype='test', fnum=fnum)

    # MANY

    def _get_centroids(self):
        cvh = CVHelper()
        cvh.set_exp('any')
        cvh.set_quality(self.quality)
        # cvh.set_quality([1])
        return get_centroids(cvh._load_results(ftype='train'))

    def _cv_many(self):
        """
        Get the maximum number of frames to require from any of the
        loaded cutout sequences. Since they have different number of
        cutout frames, it must be the lowest number of frames in a
        sequence that should be used as the largest number of frames
        to require to have a signal, consecutive or not, in order to
        label the sequence an event.

        This number N is used to create a matrix whose rows contain
        the result of one is the recorded number of signals in a frame

        """

        # List of centroids for the best accuracy
        # Get location of the best parameters (sigma, tau) for each fold
        centroids = self._get_centroids()

        # Now we have the parameter-space indices of the highest accuracy in each fold
        # Given these, we can go through each train fold again, but this time
        # saving a vector of the same length as the number of cutout frames in each sequence
        # and in the same order holding the number of signals found in each cutout frame.
        # The same goes for the test data, again with the same sigma and tau.

        for n, data in enumerate(self._folds):

            # Fold number
            fnum = n + 1

            sigma_ix, tau_ix = centroids[n]
            params = {
                'sigma': self.xp.sigmas[sigma_ix],
                'tau': self.xp.taus[tau_ix]
            }

            # Unpack training and test fold (f: features; l: labels)
            (train_f, train_l), (test_f, test_l) = data

            # Get list of signal vectors for each fold type
            signals_train = self._signals_many(train_f, **params)
            signals_test = self._signals_many(test_f, **params)

            # Save the results
            self._save_fold_results_signals(signals_train, ftype='train', fnum=fnum)
            self._save_fold_results_signals(signals_test, ftype='test', fnum=fnum)


    # Grid search
    # -----------

    # Cube of Predictions (CoP)

    # def _cube_of_predictions(self, features):
    #     """Store result of training with chosen experiment"""
    #     # self._cop = getattr(self, '_cop_' + self.xp.name)()
    #     return getattr(self, '_cop_' + self.xp.name)(features)

    # ANY

    def _cop_any(self, features):
        """
        Returns
        -------
        cop : 3D numpy.ndarray
            Predictions given each parameter combination (rows, cols)
            and for each cutout frame in the cutout sequence (depthwise).

        """

        # Shape of the cube
        shape = (self.xp.N_sigmas, self.xp.N_taus, features.shape[0])

        cop = np.zeros(shape).astype(bool)
        for cs_ix, cs in enumerate(features):

            # Set the unique file name for the predictions
            cs.set_fname_gsp(self._fname_prediction)

            # If there is already a result, and it was made after the CV began,
            # re-use it instead of re-calculating it.
            if cs.has_gs_prediction and cs.gs_prediction_time > self._cv_time:
                cop[:, :, cs_ix] = cs.gs_prediction

            else:
                cs.load()
                cs.set_quality(self.quality)
                cs.calibrate()
                cs.calculate_residuals()
                cop[:, :, cs_ix] = cs.gridsearch(**self.xp.gskw)
                cs.cleanup()

        return cop

    # BASELINE

    def _cop_baseline(self, features):
        """Return CoP for baseline experiment"""
        shape = (self.xp.N_sigmas, self.xp.N_taus, features.shape[0])
        cop = np.zeros(shape).astype(bool)
        for cs_ix, cs in enumerate(features):
            cs.set_fname_gsp(self._fname_prediction)
            if cs.has_gs_prediction and cs.gs_prediction_time > self._cv_time:
                cop[:, :, cs_ix] = cs.gs_prediction
            else:
                cop[:, :, cs_ix] = cs.gridsearch(**self.xp.gskw)
        return cop

    # MANY

    def _signals_many(self, features, **params):
        """Return list of signal vectors for experiment MANY"""
        # from pprint import pprint
        # pprint(params)
        # raise SystemExit
        signals = []
        for k, cs in enumerate(features):
            cs.set_fname_gsp(self._fname_signals(**params))
            if cs.has_gs_prediction and cs.gs_prediction_time > self._cv_time:
                print  '*' * 10 + ' Using previously saved data ' + '*' * 10
                signals.append(cs.gs_prediction.flatten())
            else:
                cs.load()
                cs.set_quality(self.quality)
                cs.calibrate()
                cs.calculate_residuals()
                signals.append(cs.experiment_many(**params))
                cs.cleanup()
        return signals

    # Matrix of accuracies

    def _moa(self, cop, labels):
        """Store training accuracies for experiment ANY.
        A.k.a.: matrix of accuracies (moa)."""
        return (
            # Broadcast along the parameter axes
            (cop == labels[None, None, :])
            # Sum along the samples (CutoutSequences)
            .sum(axis=2).astype(float)
            # Divide by the number of training samples
            / labels.size)

    # Fold I/O
    # --------

    @property
    def _opath(self):
        # Hard-coded path for now ...
        # Would be nice, if the cutout size
        # could be stored in parameters.yaml.
        opath = os.path.join(env.paths.get('cutouts'), '101x101', 'results')
        if not os.path.exists(opath):
            os.makedirs(opath)
        return opath

    @property
    def _fn_fstr(self):
        return 'moa_E-{}_Q-{}_{}_CV-{}-{}_BG-{}.csv'

    def _fn_kw(self, ftype, fnum):
        """Return current state (self._fold)"""
        return (
            self.xp.name,
            self._qstr(self.quality),
            ftype,
            self.N_folds,
            fnum,
            self._cs.get('bg_model'),
        )

    @property
    def _fn_fstr_many_signals(self):
        return 'signals_E-{}_Q-{}_{}_CV-{}-{}_BG-{}.txt'

    @property
    def _fn_fstr_many_moa(self):
        return 'moa_E-{}_Q-{}_{}_CV-{}_BG-{}.csv'

    def _fn_kw_many_moa(self, ftype):
        """Return current state (self._fold)"""
        return (
            self.xp.name,
            self._qstr(self.quality),
            ftype,
            self.N_folds,
            self._cs.get('bg_model'),
        )

    # FILENAME COMBINATIONS

    def _filenames(self, ftype):
        import itertools as it
        args = (
            [self.xp.name],
            [self._qstr(self.quality)],
            [ftype],
            [self.N_folds],
            range(1, self.N_folds + 1),
            [self._cs.get('bg_model')],
        )
        for t in it.product(*args, repeat=1):
            ifname = os.path.join(self._opath, self._fn_fstr.format(*t))
            print ifname
            yield ifname

    def _filenames_signals(self, ftype):
        import itertools as it
        args = (
            [self.xp.name],
            [self._qstr(self.quality)],
            [ftype],
            [self.N_folds],
            range(1, self.N_folds + 1),
            [self._cs.get('bg_model')],
        )
        for t in it.product(*args, repeat=1):
            ifname = os.path.join(self._opath, self._fn_fstr_many_signals.format(*t))
            print ifname
            yield ifname

    # SAVE

    def _save_fold_results_moa(self, moa, ftype, fnum):
        fname = self._fn_fstr.format(*self._fn_kw(ftype=ftype, fnum=fnum))
        ofname = os.path.join(self._opath, fname)
        print 'Saving fold results to:', fname
        df = pd.DataFrame(data=moa, index=self.xp.sigmas, columns=self.xp.taus)
        df.to_csv(ofname, index=True, header=True)
        return ofname

    def _save_fold_results_signals(self, signals, ftype, fnum):
        fname = self._fn_fstr_many_signals.format(*self._fn_kw(ftype=ftype, fnum=fnum))
        ofname = os.path.join(self._opath, fname)
        print 'Saving fold results to:', fname
        with open(ofname, 'w+') as fsock:
            for signal_vector in signals:
                fsock.write(','.join(signal_vector.astype(str)) + '\n')
        return ofname

    # LOAD

    def _load_results(self, ftype):
        """
        Based on given parameters and chosen experiment, load all
        combinations of filenames possible with the given parameters.

        """
        #     getattr(self, '_load_results_' + self.xp.name)()

        # def _load_results_any(self, ftype):
        #     """Load results from experiment ANY"""

        # Cube of accuracies
        N_folds = self.N_folds
        coa = np.zeros((self.xp.N_sigmas, self.xp.N_taus, N_folds))
        for fold_ix, ifname in enumerate(self._filenames(ftype)):
            df = pd.read_csv(ifname, index_col=[0])
            coa[:, :, fold_ix] = df.values
        return coa

    def _load_results_signals(self, ftype):
        """
        Based on given parameters and chosen experiment,
        load signal vectors for each fold and fold type.

        Returns
        -------
        folds : list
            List with as many entries as folds.
            Each list entry is another list; this one containing signal
            vectors; one for each cutout sequence in the given fold.

        """
        folds = []
        for ifname in self._filenames_signals(ftype):
            with open(ifname, 'r') as fsock:
                folds.append(lines2intarrays(fsock.readlines()))
        return folds

    def _load_prediction_cubes(self, ftype=None):
        """
        Return prediction cube for each fold and given the fold type `ftype`.
        """

        if ftype is None:
            return

        # Lists of cubes and labels (vectors) for those cubes
        cops = []
        labels = []

        # Unpack training and test fold (f: features; l: labels)
        for (train_f, train_l), (test_f, test_l) in self._folds:

            # Use only the part of the data that is asked for.
            if ftype == 'test':
                features = test_f
                labels.append(test_l)

            elif ftype == 'train':
                features = train_f
                labels.append(train_l)

            # Create new cube-of-predictions container
            shape = (self.xp.N_sigmas, self.xp.N_taus, features.shape[0])
            cop = np.zeros(shape).astype(bool)

            for cs_ix, cs in enumerate(features):
                cs.set_fname_gsp(self._fname_prediction)
                cop[:, :, cs_ix] = cs.gs_prediction

            cops.append(cop)

        return cops, labels

    # Visualisation
    # -------------

    def plot(self):
        """Plot results for chosen experiment"""
        getattr(self, '_plot_results_' + self.xp.name)()

    def _plot_results_any(self): self._plot_results_accuracies()
    def _plot_results_any2(self): self._plot_results_accuracies()

    def _plot_results_blr(self): self._plot_results_accuracies()
    def _plot_results_bla(self): self._plot_results_accuracies()
    def _plot_results_bln(self): self._plot_results_accuracies()

    def _plot_results_blr2(self): self._plot_results_accuracies()
    def _plot_results_bla2(self): self._plot_results_accuracies()
    def _plot_results_bln2(self): self._plot_results_accuracies()

    def _plot_results_accuracies(self):
        """Plot results for experiment ANY and BASELINE experiments"""

        # from scipy import stats
        import matplotlib as mpl
        # Select backend beforehand to make it run on the image servers
        # Run anywhere, where matplotlib has already been loaded, this
        # will produce a warning that the backend has already been chosen.
        # Just ignore this.
        mpl.use('pdf')
        import matplotlib.pyplot as plt

        from mplconf import mplrc
        from mplconf import rmath

        # Load data (cube of accuracies (COA))
        coa_train = self._load_results(ftype='train')
        coa_test = self._load_results(ftype='test')

        # Get CV information
        N_folds = self.N_folds

        # Get experiment information
        qstr = self._qstr(self.quality)

        # Extract parameter space

        """
        NOTE
        ----
        All though these values are extracted from the parameter file,
        the cross validation may not have been run for these values of
        sigma and tau.
            Alternatively, the code should take into account what range
        of hyper-parameters were used (saved along with the results).
        """

        sigmas = self.xp.sigmas
        taus = self.xp.taus

        N_sigmas = sigmas.size
        N_taus = taus.size

        # Prepare the data for plotting
        # -----------------------------

        # Plot of every Matrix of Accuracy (MOA) for each training fold
        # Check: This is what is loaded directly from the files on disk.

        # On top of this plot, I need the locations of the best
        # accuracies all together.
        # Matrices of maximum accuracies (momas)

        """True => accuracy for parameter-pair location is max."""

        # Checked function `get_matrices_of_maximum_accuracy()` through
        # on 2014-06-14, and it looked ok.
        momas_train = get_matrices_of_maximum_accuracy(coa_train)
        momas_test = get_matrices_of_maximum_accuracy(coa_test)

        """
        Calculate the locations of the assumed best accuracies
        This is the CENTRE-OF-MASS (COM)/CENTROID for each set of
        entries with the maximum value.
        """

        # List of centroids for the best accuracy

        centroids_train = get_centroids(coa_train)
        centroids_test = get_centroids(coa_test)

        # Use slices to create matrices with NaNs
        # everywhere except at the centroid index.

        shape = (N_sigmas, N_taus)
        coms_train = get_centroid_stack(centroids_train, shape, fill=np.nan)
        coms_test = get_centroid_stack(centroids_test, shape, fill=np.nan)

        # # Density estimation of the highest-accuracy locations
        # kernel_train = stats.gaussian_kde(coms_train.sum(axis=2))

        # # Density estimation of the highest-acuracies for the sum of the test accuracies.
        # kernel_test = stats.gaussian_kde(coa_test.sum(axis=2))

        # Plot settings
        # -------------

        # mplrc('publish_digital')
        mplrc('publish_printed')

        def imshow_com(img, ax=None):
            """Plot centre of mass image"""
            if ax is None: ax = plt.gca()
            ix = (img == 1)
            com = np.zeros((N_sigmas, N_taus, 4))
            com[ix, :] = np.array([1., .0, .0, 1.])
            com[~ix, :] = np.array([0., .0, .0, .0])
            ax.imshow(com, **comskw)

        # Extent
        # With the bottom is sigmas.min(), imshow must have origin='lower'
        tmin, tmax = taus.min(), taus.max()
        dtau = (taus[1] - taus[0]) / 2.
        smin, smax = sigmas.min(), sigmas.max()
        dsig = (sigmas[1] - sigmas[0]) / 2.
        extent = [(tmin - dtau), (tmax + dtau), (smin - dsig), (smax + dsig)]

        # GENERIC image settings
        imkw = dict(origin='lower', aspect='auto', interpolation='nearest')
        imkw.update(extent=extent)

        # Image settings for matrices of accuracies
        moaskw = imkw.copy()
        moaskw.update(cmap=mpl.cm.coolwarm)
        # Accuracies can only be between 0 and 100 percent
        # vmin = 0
        # vmax = 1
        vmin = .4
        vmax = .6
        moaskw.update(vmin=vmin, vmax=vmax)

        # Image settings for matrices of maximum-accuracy indices
        momaskw = imkw.copy()
        momaskw.update(cmap=mpl.cm.binary)

        # Image settings for matrices of centre-of-mass index
        comskw = imkw.copy()
        comskw.update(cmap=None)

        # Figure settings (for plt.subplots)
        ncols = 5
        nrows = N_folds / ncols
        figh = [3.3, 5.7]
        figsize = (13.5, figh[nrows - 1])
        fkw = dict(sharex=True, sharey=True, figsize=figsize)
        # fkw = dict(figsize=(15, 6.5))

        # Figure-adjustement settings
        figb = [.15, .15]
        adjustkw = dict(left=.04, bottom=figb[nrows - 1], top=.82, right=.95)
        adjustkw.update(wspace=.10, hspace=.10)

        # Colorbar settings
        cbkw = dict(extend='neither', drawedges=False)
        # cbkw.update(ticks=np.linspace(.0, 1., 3))
        cbkw.update(ticks=np.linspace(vmin, vmax, 3))
        cbkw.update(orientation='vertical')
        # rect = [left, bottom, width, height]
        # crect = [0.96, 0.52, 0.01, 0.43]
        rect_b = adjustkw['bottom']
        rect_h = adjustkw['top'] - adjustkw['bottom']
        crect = [0.96, rect_b, 0.01, rect_h]

        # Default title for all plots
        qtit = rmath('Quality combination: {}'.format(qstr))
        titkw = dict(fontsize=18)

        # Labels
        lblkw = dict(fontsize=15)

        # Accuracies
        # ----------

        """TRAIN: Plot accuracies for each training fold"""

        # Create the figure and axes
        fig, axes = plt.subplots(nrows, ncols, **fkw)

        # axes.flat[ncols / 2].set_title(qtit, **titkw)
        fig.suptitle(qtit, **titkw)

        # Plot data on eaxh axis
        for i, ax in enumerate(axes.flat):
            im = ax.imshow(coa_train[:, :, i], **moaskw)
            ax2 = ax.twiny()
            ax2.set_xticks([])
            ax2.set_xlabel(rmath('Fold {:d}'.format(i + 1)), **lblkw)

        scl_kw = dict(l=r'$\sigma$', b=r'$\tau$', bticks=[.0, .5, 1.])
        scl_kw.update(lblkw=lblkw)
        set_common_labels(axes, ncols, **scl_kw)
        ax.set_xticks([.0, .5, 1.])

        # Plot colorbar to the right of the accuracies
        cax = plt.axes(crect)
        fig.colorbar(mappable=im, cax=cax, **cbkw)

        # fig.tight_layout()
        fig.subplots_adjust(**adjustkw)
        fname = 'moa_E-{}_Q-{}_CV-{}_train_folds.pdf'.format(
            self.xp.name, qstr, N_folds)
        ofname = os.path.join(self._opath, fname)
        plt.savefig(ofname)
        plt.close(fig)


        """TRAIN: Plot best-accuracy positions with centroids"""

        # Create the figure and axes
        fig, axes = plt.subplots(nrows, ncols, **fkw)

        # axes.flat[ncols / 2].set_title(qtit, **titkw)
        fig.suptitle(qtit, **titkw)

        # Plot data on eaxh axis
        for i, ax in enumerate(axes.flat):
            # Show the locations of all entries having the maximum accuracy
            ax.imshow(momas_train[:, :, i], **momaskw)
            ax2 = ax.twiny()
            ax2.set_xticks([])
            ax2.set_xlabel(rmath('Fold {:d}'.format(i + 1)), **lblkw)

            # Plot a red color for the location of the centre-of-mass
            # of the maximum accuracy indices.
            imshow_com(img=coms_train[:, :, i], ax=ax)

        set_common_labels(axes, ncols, **scl_kw)
        ax.set_xticks([.0, .5, 1.])

        # fig.tight_layout()
        fig.subplots_adjust(**adjustkw)
        fname = 'moa_E-{}_Q-{}_CV-{}_train_best.pdf'.format(
            self.xp.name, qstr, N_folds)
        ofname = os.path.join(self._opath, fname)
        plt.savefig(ofname)
        plt.close(fig)

        """TEST: Plot accuracies for each test fold"""

        # Create the figure and axes
        fig, axes = plt.subplots(nrows, ncols, **fkw)

        # axes.flat[ncols / 2].set_title(qtit, **titkw)
        fig.suptitle(qtit, **titkw)

        # Plot data on eaxh axis
        for i, ax in enumerate(axes.flat):
            im = ax.imshow(coa_test[:, :, i], **moaskw)
            ax2 = ax.twiny()
            ax2.set_xticks([])
            ax2.set_xlabel(rmath('Fold {:d}'.format(i + 1)), **lblkw)

        # Set the x label on the bottom-most axes only
        for ax in axes.flat[-ncols:]:
            ax.set_xlabel(r'$\tau$')
            ax.set_xticks([.0, .5, 1.])

        # Set the y label on the left-most axes only
        for ax in axes.flat[::ncols]:
            ax.set_ylabel(r'$\sigma$')

        # Plot colorbar to the right of the accuracies
        cax = plt.axes(crect)
        fig.colorbar(mappable=im, cax=cax, **cbkw)

        # fig.tight_layout()
        fig.subplots_adjust(**adjustkw)
        fname = 'moa_E-{}_Q-{}_CV-{}_test_folds.pdf'.format(
            self.xp.name, qstr, N_folds)
        ofname = os.path.join(self._opath, fname)
        plt.savefig(ofname)
        plt.close(fig)

        """TEST: Plot best-accuracy positions with centroids"""

        # Create the figure and axes
        fig, axes = plt.subplots(nrows, ncols, **fkw)

        # axes.flat[ncols / 2].set_title(qtit, **titkw)
        fig.suptitle(qtit, **titkw)

        # Plot data on eaxh axis
        for i, ax in enumerate(axes.flat):
            # Show the locations of all entries having the maximum accuracy
            ax.imshow(momas_test[:, :, i], **momaskw)
            ax2 = ax.twiny()
            ax2.set_xticks([])
            ax2.set_xlabel(rmath('Fold {:d}'.format(i + 1)), **lblkw)

            # Plot a red color for the location of the centre-of-mass
            # of the maximum accuracy indices.
            imshow_com(img=coms_test[:, :, i], ax=ax)

        set_common_labels(axes, ncols, **scl_kw)
        ax.set_xticks([.0, .5, 1.])

        # fig.tight_layout()
        fig.subplots_adjust(**adjustkw)
        fname = 'moa_E-{}_Q-{}_CV-{}_test_best.pdf'.format(
            self.xp.name, qstr, N_folds)
        ofname = os.path.join(self._opath, fname)
        plt.savefig(ofname)
        plt.close(fig)

        """EXIT"""
        """EXIT"""
        """EXIT"""
        return

        """
        Plot distribution of the locations of
        the best accuracies in the training folds.
        """

        # fkw = dict(sharex=True, sharey=True, figsize=(15, 8))
        # fig, (ax1, ax2) = plt.subplots(1, 2, **fkw)

        # X, Y = np.mgrid[0:1:N_sigmas * 1j, 0:1:N_taus * 1j]
        # positions = np.vstack([X.ravel(), Y.ravel()])
        # Z = np.reshape(kernel_train(positions).T, X.shape)
        # ax1.imshow(np.rot90(Z), **moaskw)

        # Plot colorbar above the accuracies...
        # ax = fig.add_axes([.25, 1.2, .5, .0])
        # fig.colorbar(mappable=im, ax=ax, **cbkw)

        # fig.tight_layout()
        # fname = 'kde_coms+coa_test.sum.pdf'
        # ofname = os.path.join(self._opath, fname)
        # plt.savefig(ofname)

        """
        Plot distribution of
        the sum of the test-fold accuracy matrices.
        """

        """
        Plot the surface of
        average accuracy for the test folds
        along with standard-deviation surfaces.
        """

    def _plot_results_many(self):
        """Plot results for experiment MANY"""

        import matplotlib as mpl
        # Select backend beforehand to make it run on the image servers
        mpl.use('pdf')
        import matplotlib.pyplot as plt

        from mplconf import mplrc
        from mplconf import rmath

        mplrc('publish_digital')

        # Load data

        # Train
        fname_train = self._fn_fstr_many_moa.format(*self._fn_kw_many_moa(ftype='train'))
        ifname_train = os.path.join(self._opath, fname_train)
        moa_train = pd.read_csv(ifname_train, index_col=[0])

        # Number of frames required at the chosen maximum
        train_acc_max_ix = np.argmax(moa_train.values, axis=1)

        # Test
        fname_test = self._fn_fstr_many_moa.format(*self._fn_kw_many_moa(ftype='test'))
        ifname_test = os.path.join(self._opath, fname_test)
        moa_test = pd.read_csv(ifname_test, index_col=[0])

        # Get the vector with numbers of required frames
        # N_frames = moa_train.columns.astype(int)
        N_frames = np.arange(0, moa_train.shape[1])

        # Get the maximum number of frames that can
        # be required when using this quality
        cs_frame_count = np.zeros(self._css.size).astype(int)
        for cs_ix, cs in enumerate(self._css):
            cs.set_quality(self.quality)
            # cs.calibrate()
            cs_frame_count[cs_ix] = len(cs)
        N_max_frames = np.min(cs_frame_count)

        # Get CV information
        N_folds = self.N_folds
        # Get experiment information
        qstr = self._qstr(self.quality)

        # Extract chosen parameters
        centroids = self._get_centroids()

        # Plot settings
        # -------------

        colors = mpl.rcParams.get('axes.color_cycle')

        # pkw = dict(lw=3, c=colors[0])
        inkw = dict(ls='none', marker='.', ms=8, mec='none', c=colors[4])
        # Make the point visible outside the axes extent,
        # and put it on top of everythin else in the axes instance
        inkw.update(clip_on=False, zorder=100)
        # bboxkw = dict(facecolor='#efefef', edgecolor='#cccccc', pad=10.0)
        # bboxkw = dict()
        bboxkw = dict(facecolor='w', alpha=.8)
        tkw = dict(fontsize=12, ha='left', va='bottom', bbox=bboxkw)
        fstr1 = r'\sigma = {:.2f}'
        fstr2 = r'\tau = {:.2f}'
        fstr3 = r'Train max (\nu = {: >2d})'
        fstr4 = r'Largest \nu for given quality combo (\nu = {: >2d})'
        # fstr = r'\sigma = {:.2f}' + '\n' + r'\tau = {:.2f}'
        s4 = fstr4.format(N_max_frames)

        # Labels
        xlabel = rmath(r'\nu / Minimum required number of frames with a signal')

        # Accuracies
        # ----------

        # Create the figure and axes
        fkw = dict(sharex=True, sharey=True, figsize=(13., 8.))
        fig, axes = plt.subplots(N_folds, 1, **fkw)

        # Plot data on eaxh axis
        insert_data = []
        for fold_ix, ax in enumerate(axes.flat):

            if fold_ix == 0:
                ax.set_title(rmath('Quality combination: {}'.format(qstr)), fontsize=18)

            # Add the accuracy lines

            # ax.plot(N_frames, moa_train.values[fold_ix, :], **pkw)
            ax.plot(N_frames, moa_train.values[fold_ix, :], label=rmath('Train'))
            ax.plot(N_frames, moa_test.values[fold_ix, :], label=rmath('Test'))
            N_frames_best = train_acc_max_ix[fold_ix]
            s3 = fstr3.format(N_frames_best)
            ax.axvline(x=N_frames_best, c='k', label=rmath(s3))
            ax.axvline(x=N_max_frames, c='r', label=rmath(s4))
            # print ax.bbox

            # Display the values of \sigma and \tau

            fnum = fold_ix + 1
            sigma_ix, tau_ix = centroids[fold_ix]
            sigma, tau = self.xp.sigmas[sigma_ix], self.xp.taus[tau_ix]

            # Save handles and data until after tightening the figure layout
            insert_data.append((ax, sigma, tau))

            # Write data out immediately
            s1 = rmath(fstr1.format(sigma))
            s2 = rmath(fstr2.format(tau))
            # s = '{}\n\{}'.format(s1, s2)
            # s = rmath(fstr.format(sigma, tau))
            # y = ax.get_position().ymin + .01
            ax.text(.925, .35, s1, transform=ax.transAxes, **tkw)
            ax.text(.925, .15, s2, transform=ax.transAxes, **tkw)

            ax.legend(ncol=3)
            ax2 = ax.twinx()
            ax.set_ylabel(rmath('Accuracy'))
            ax2.set_ylabel(rmath('Fold {}'.format(fnum)))

            # ax.set_yticks([.0, .5, 1.])
            ax2.set_yticks([])

        ax.set_xlabel(xlabel)
        ax.set_xlim(np.min(N_frames), np.max(N_frames))
        # ax.set_ylim(.0, 1.)

        fig.tight_layout()

        # Add insert axes (after the change in fig.transFigure from calling fig.tight_layout())
        fig_w, fig_h = fig.get_size_inches()
        fig_w2h = float(fig_w) / fig_h
        sz = .035
        l, w, h = .85, sz, sz * fig_w2h
        for ax, sigma, tau in insert_data:
            b = ax.get_position().ymin + .015
            axin = fig.add_axes([l, b, w, h])
            axin.plot([tau], [sigma], **inkw)
            # axin.set_ylabel(rmath(r'\sigma'))
            # axin.set_xlabel(rmath(r'\tau'))
            axin.set_ylim(*self.xp.sigmas[[0, -1]])
            axin.set_xlim(*self.xp.taus[[0, -1]])
            axin.set_xticks([])
            axin.set_yticks([])

        fname = 'moa_E-{}_Q-{}_CV-{}.pdf'.format(self.xp.name, qstr, N_folds)
        ofname = os.path.join(self._opath, fname)
        plt.savefig(ofname)
        plt.close(fig)

    # Analysis
    # --------

    def analyse(self):
        """Analyse results for chosen experiment"""
        getattr(self, '_analyse_' + self.xp.name)()

    def _analyse_any(self): self._analyse_simple()
    def _analyse_any2(self): self._analyse_simple()

    def _analyse_blr(self): self._analyse_simple() # baseline()
    def _analyse_bla(self): self._analyse_simple() # baseline()
    def _analyse_bln(self): self._analyse_simple() # baseline()

    def _analyse_blr2(self): self._analyse_simple() # baseline()
    def _analyse_bla2(self): self._analyse_simple() # baseline()
    def _analyse_bln2(self): self._analyse_simple() # baseline()

    def _analyse_simple(self):
        """Analyse results any experiment"""

        # print self._fname_prediction
        # raise SystemExit

        # Load cubes of accuracies
        # (one matrix-of-accuracy for each fold and fold type)
        coa_train = self._load_results(ftype='train')

        # List of centroids for the best accuracy

        # Get location of the best parameters for each fold
        centroids = get_centroids(coa_train)

        # Use these indices to extract the predictions for each test fold
        # ultimately using these to get, not just accuracies, but the entire
        # table of confusion based on the confusion matrix (which for binary
        # classification is the same).
        # Get list of prediction cubes, one for each fold
        cops, labels = self._load_prediction_cubes(ftype='test')

        # print labels
        # raise SystemExit

        # Get a list of vectors containing the predictions
        predictions = []
        for fold_ix, (sigma_ix, tau_ix) in enumerate(centroids):
            predictions.append(cops[fold_ix][sigma_ix, tau_ix, :])

        display_table_of_confusion(predictions, labels)

    def _analyse_many(self):

        folds_train = self._load_results_signals(ftype='train')
        folds_test = self._load_results_signals(ftype='test')

        N_max_frames = np.min([len(cs) for cs in self._css])
        print 'Maximum number of needed frames required (using all quality combinations):', N_max_frames

        # Include zero (always yes) and the max number of frames
        N_frames = np.arange(0, N_max_frames + 1)

        # Cube of predictions (USED TO CALCULATE THE TABLE OF CONFUSION)
        #  Depth-wise: fold index
        #  Row-wise: cutout sequence
        #  Column-wise: Prediction given required number of frames with signals
        #  Number of columns is number of frames to require.
        N_css = len(self._css)
        N_test = N_css // self.N_folds
        N_train = N_css - N_test
        print 'N_train', N_train
        print 'N_test', N_test
        print 'N_train + N_test', N_train + N_test
        cop_train = np.zeros((N_train, N_frames.size, self.N_folds)).astype(bool)
        cop_test = np.zeros((N_test, N_frames.size, self.N_folds)).astype(bool)

        # Matrix of accuracies (Saved for later and used for plotting)
        # One row for each fold
        # Number of columns is number of frames to require.
        moa_train = np.zeros((self.N_folds, N_frames.size))
        moa_test = np.zeros((self.N_folds, N_frames.size))

        # Experiment
        labels = []
        z = zip(folds_train, folds_test, self._folds)
        for fold_ix, (train, test, data) in enumerate(z):

            # Unpack training and test fold (f: features; l: labels)
            (train_f, train_l), (test_f, test_l) = data

            # Add the test labels for later
            labels.append(test_l)

            # Create a matrix of predictions
            mop_train = np.zeros((N_train, N_frames.size)).astype(bool)
            mop_test = np.zeros((N_test, N_frames.size)).astype(bool)

            # TRAIN

            # Train (continuing where ANY left off)
            for train_ix, signal_vector in enumerate(train):
                N_signals = (signal_vector > 0).sum()
                # Vector of predictions
                # If number of frames containing (any number of) signals
                # is at least as large as the number of required frames,
                # predict that it is an event.
                mop_train[train_ix, :] = (N_frames <= N_signals)

            # Save the matrix of predictions in the cube of predictions
            cop_train[:, :, fold_ix] = mop_train

            # Vector of accuracies for given fold
            # Add it as a row in the matrix of accuracies.
            moa_train[fold_ix, :] = (
                # Broadcast condition column wise
                (mop_train == train_l[:, None])
                # Sum over training samples
                .sum(axis=0).astype(float)
                # Divide by total number of samples
                / N_train)

            # TEST

            # Train (continuing where ANY left off)
            for test_ix, signal_vector in enumerate(test):
                N_signals = (signal_vector > 0).sum()
                # Vector of predictions
                # If number of frames containing (any number of) signals
                # is at least as large as the number of required frames,
                # predict that it is an event.
                mop_test[test_ix, :] = (N_frames <= N_signals)

            # Save the matrix of predictions in the cube of predictions
            cop_test[:, :, fold_ix] = mop_test

            # Vector of accuracies for given fold
            # Add it as a row in the matrix of accuracies.
            moa_test[fold_ix, :] = (
                # Broadcast condition column wise
                (mop_test == test_l[:, None])
                # Sum over test samples
                .sum(axis=0).astype(float)
                # Divide by total number of samples
                / N_test)

        # Save the accuracies
        # -------------------

        folds_ix = np.arange(self.N_folds) + 1

        # Train
        fname_train = self._fn_fstr_many_moa.format(*self._fn_kw_many_moa(ftype='train'))
        ofname_train = os.path.join(self._opath, fname_train)
        print ofname_train
        df = pd.DataFrame(data=moa_train, index=folds_ix, columns=N_frames)
        df.to_csv(ofname_train, index=True, header=True)

        # Test
        fname_test = self._fn_fstr_many_moa.format(*self._fn_kw_many_moa(ftype='test'))
        ofname_test = os.path.join(self._opath, fname_test)
        print ofname_test
        df = pd.DataFrame(data=moa_test, index=folds_ix, columns=N_frames)
        df.to_csv(ofname_test, index=True, header=True)

        # Report some numbers
        # -------------------

        # Number of frames giving highest accuracy for each fold

        if 0:
            # Get the maximum accuracy for each fold.
            train_acc_max = moa_train.max(axis=1)

            # Given the maximum accuracy in each fold,
            train_acc_max_ix = np.argmax(moa_train == train_acc_max[:, None], axis=1)

        else:
            # This yields the same result
            train_acc_max_ix = np.argmax(moa_train, axis=1)

        # Print out the maximum accuracy of the training and the number of frames required to get it.
        print 'TRAINING:'
        print '\n'.join([
            ' Fold {:d}: Best acc.: {7.5f}; N frames: {:d}'.format(
                fix, moa_train[fix, Nix], Nix) for fix, Nix in enumerate(
                    train_acc_max_ix)])

        # I need a confusion matrix for each fold
        # For that I need a prediction vector for each test fold
        # To make it, I need the best choice of number of frames to
        # require for each fold and then get the predictions given this choice.

        # The best-number-of-frames index is along the columns in the cube of predictions,
        # (So this is a different form when I obtain the predictions in the previous experiment.)
        # For each fold (fold_ix), I take the best-accuracy index (N_frames_ux) at which the number of required signal frames have the highest accuracy
        # and grab the prediction vector cop_test[:, N_frames_ix, fold_ix] (i.e. along the rows of the cop)
        # and calculate the confusion matrix given the labels of the test set

        # Prepare input for `display_table_of_confusion`.
        predictions = []
        for fold_ix in range(self.N_folds):

            N_frames_ix = train_acc_max_ix[fold_ix]
            predictions.append(cop_test[:, N_frames_ix, fold_ix])

        display_table_of_confusion(predictions, labels)

    # -- END class CVHelper --


def display_table_of_confusion(predictions, labels):

    # Get cube of confusion matrices
    classes = np.array([True, False])
    coc = confusion_cube(predictions, labels, classes)
    print 'coc.shape =', coc.shape
    # print 'coc[:, :, 0]:'
    # print coc[:, :, 0]

    # print labels[0]
    # print predictions[0]
    # raise SystemExit

    mean = coc.mean(axis=2)
    std = coc.std(axis=2)
    err = std / np.sqrt(coc.shape[2])

    # count_vector contains total number of entries in each fold
    count_vector = coc.sum(axis=0).sum(axis=0).astype(float)

    # Divide by total number of entries for given fold to get accuracy
    #   Assuming here that the class of interest is the one indexed
    #   first along the rows and columns of the confusion matrix,
    #   i.e. the class SN.
    acc_vector = (coc[0, 0, :] + coc[1, 1, :]) / count_vector
    acc_mean = acc_vector.mean()
    acc_std = acc_vector.std()
    acc_err = acc_std / np.sqrt(coc.shape[2])

    # Print all results
    print '\nClass order in confusion matrix:', list(classes)

    print '\nMean table of confusion:'
    print mean
    print '\nStandard deviation:'
    print std
    print '\nStandard error:'
    print err

    print '\nMean accuracy (TP + TN):', acc_mean,
    print '+-', acc_err, '(std = {})'.format(acc_std)


def max_indices(arr):
    """Get indices of the maximum"""
    return (arr == arr.max()).nonzero()


def max_locs(arr):
    """
    Get boolean array where True entries
    mark location of maximum value
    """
    return (arr == arr.max())


def get_matrices_of_maximum_accuracy(cube):
    """
    Get stack of boolean matrices where
    True marks location of maximum value
    for each layer in the stack.
    """
    # momas = np.zeros_like(cube).astype(bool)
    # for i in range(cube.shape[2]):
    #     momas[:, :, i] = max_locs(cube[:, :, i])
    # return momas
    return (cube == cube.max(axis=0).max(axis=0)[None, None, :])


def centroid(indices_list):
    """Return indices nearest to centre-of-mass of given indices"""
    return [np.round(np.mean(indices)) for indices in indices_list]


def get_centroids(cube):
    """Return list of centroids depth-wise for a data cube"""
    return [centroid(max_indices(cube[:, :, fold_ix]))
            for fold_ix in range(cube.shape[2])]


def get_centroid_stack(centroids, shape, fill=np.nan):
    """
    Create data cube with ones at the location
    of the centroids and `fill` elsewhere.
    """
    coms = np.zeros(shape + (len(centroids),))
    for i, (rix, cix) in enumerate(centroids):
        coms[rix, cix, i] = 1
    coms[coms == 0] = fill
    return coms


def set_common_labels(axes, ncols,
    l=None, b=None,
    lticks=None, bticks=None,
    lblkw={}):

    if b is not None:
        # Set the x label on the bottom-most axes
        for ax in axes.flat[-ncols:]:
            ax.set_xlabel(b, **lblkw)

    if l is not None:
        # Set the y label on the left-most axes
        for ax in axes.flat[::ncols]:
            ax.set_ylabel(l, **lblkw)


def load_parameters(ifname=env.files.get('params')):
    """Load parameters for loading the data and running an experiment"""
    with open(ifname, 'r') as fsock:
        return yaml.load(fsock.read())


def dict2numpy(d):
    """
    Walk through given dictionary recursively.

    If any sub-dictionary has keys that correspond to argument names
    for NunmPy functions such as

        numpy.linspace(start, stop, [num])  # here: start, stop, num
        numpy.arange([start,] stop[, step]) # here: start, stop, step

    the dictionary is replaced by the corresponding construct, i.e. either of

        numpy.linspace(**d)
        numpy.arange(**d)

    All three arguments must be specified in the dictionary for the parser to
    manipulate them.

    Parameters
    ----------
    d : dict
        A dictionary whose keys are to be parsed recursively.

    Returns
    -------
    d : dict
        A dictionary with any sub-dictionaries parsed recursively.
    L : numpy.linspace(**d)
    A : numpy.arange(**d)

    """

    if 3 == len(d.keys()):

        keys = ('start', 'stop',)

        if 3 == sum([d.has_key(k) for k in keys + ('num',)]):
            return np.linspace(**d)

        elif 3 == sum([d.has_key(k) for k in keys + ('step',)]):
            return np.arange(**d)

    for key, val in d.items():
        if isinstance(val, dict):
            d[key] = dict2numpy(val)

    return d


def get_folds(features, labels, N_folds):
    """Generator for easy loop over folds in cross validation."""

    features = np.array(features)
    labels = np.array(labels)

    # This seems reasonable to assume
    assert features.shape[0] == labels.shape[0]

    # For this project, I make sure that the
    # number of entries is divisible by `N_folds`
    assert features.shape[0] % N_folds == 0

    # Define chunk sizes
    # ------------------

    # Create training and test data sets
    N_features = features.shape[0]
    # Number of items in the test set
    # N_items is the number of samples per chunk
    N_items = N_features / N_folds

    # Define slices that will cut out the parts for training and testing
    # -------------
    tests = [
        slice(i * N_items, (i + 1) * N_items, None)
        for i in range(N_folds)]
    trains = [(
        slice(0, i * N_items, None),
        slice((i + 1) * N_items, N_folds * N_items, None))
        for i in range(N_folds)]

    # Generate the iterator
    for (s, (s1, s2)) in zip(tests, trains):

        test = (features[s], labels[s])

        train_features = np.append(features[s1], features[s2])
        train_labels = np.append(labels[s1], labels[s2])
        train = (train_features, train_labels)

        yield train, test


def lines2intarrays(lines):
    list_ = []
    for line in lines:
        list_.append(np.array(line.split(',')).astype(int))
    return list_


# -----------------------------------------------------------------------------

def table_of_confusion():
    """
                  Prediction
                  ----------

                  T         F

             ---------- ----------
    A|      |          |          |
    c|  T   |    TP    |    FN    |
    t|      |          | (Type 2) |
    u|       ---------- ----------
    a|      |          |          |
    l|  F   |    FP    |    TN    |
            | (Type 1) |          |
             ---------- ----------

    """
    pass


def _confusion_matrix(predictions, labels, classes=None):
    """

    Parameters
    ----------
    classes : 1D-array
        List of classes to compare against.
        If provided, the order of the entries is preserved along
        the rows and columns of the returned confusion matrix.

    Returns
    -------
    CFM : 2D-array
        The confusion matrix of shape (N, N) for N classes.
        The confusion matrix corresponding to every possible
        class that appears in `labels`.

    """

    import itertools as it

    predictions = np.asarray([predictions]).flatten()
    labels = np.asarray([labels]).flatten()
    if classes is None:
        # np.unique also sorts the data
        classes = np.unique(labels)
    else:
        pass
        # Can not use np.unique, since this also sorts the input array.
        # Can not use Python's set() function since its elements are unordered.
        # classes = np.array(list(set(classes)))
    N = classes.size

    # The order (0, 1) of the arguments for zip() translates into
    # (row-wise, col-wise) location of the two inputs. In this case,
    # since labels are given before predictions, the labels will be seen
    # along the rows, and the predictted classes shown  along the columns.
    z = zip(labels, predictions)
    p = it.product(classes, repeat=2)
    # print list(p)
    return np.array([z.count(x) for x in p]).reshape(N, N)


def confusion_cube(predictions, labels, classes=None):
    from sklearn.metrics import confusion_matrix as cm
    return np.dstack([cm(l, p, classes) for (p, l) in zip(predictions, labels)])


# -----------------------------------------------------------------------------


def cv(exp='any', quality=None):
    """Perform N-fold cross-validated grid search for chosen experiment."""
    cvh = CVHelper()
    cvh.set_exp(exp)
    if quality is not None:
        cvh.set_quality(quality)
    cvh.cv()


def analyse(exp='any', quality=None):
    """Analyse results of N-fold cross-validated grid search for chosen experiment."""
    cvh = CVHelper()
    cvh.set_exp(exp)
    if quality is not None:
        cvh.set_quality(quality)
    cvh.analyse()


def plot(exp='any', quality=None):
    """Plot results of N-fold cross-validated grid search for chosen experiment."""
    cvh = CVHelper()
    cvh.set_exp(exp)
    if quality is not None:
        cvh.set_quality(quality)
    cvh.plot()
