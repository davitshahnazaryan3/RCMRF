"""
Postprocessing of results
"""
from pathlib import Path
import pickle
import numpy as np
import pandas as pd
import os
import json
from scipy.interpolate import interp1d
from scipy.stats import lognorm


class IDAPostprocessor:
    def __init__(self, path, export=True, export_path=None, flag3d=False):
        """
        :param path: str                    IDA results directory
        :param export: bool                 Exporting postprocessed IDA results or not
        :param flag3d: bool                 Run IDA for 3D or not
        """
        self.path = path
        self.export = export
        self.export_path = export_path
        self.flag3d = flag3d

    def export_results(self, filepath, data, filetype):
        """
        Store results in the database
        :param filepath: str                            Filepath, e.g. "directory/name"
        :param data:                                    Data to be stored
        :param filetype: str                            Filetype, e.g. npy, json, pkl, csv
        :return: None
        """
        if filetype == "npy":
            np.save(f"{filepath}.npy", data)
        elif filetype == "pkl" or filetype == "pickle":
            with open(f"{filepath}.pickle", 'wb') as handle:
                pickle.dump(data, handle)
        elif filetype == "json":
            with open(f"{filepath}.json", "w") as json_file:
                json.dump(data, json_file)
        elif filetype == "csv":
            data.to_csv(f"{filepath}.csv", index=False)

    def cscvn2(self, points):
        cs = interp1d(points[0], points[1])
        return cs

    def splinequery_IDA(self, spl, edp_range, qtile_range, edp, num_recs):

        # Getting the edp-based im values
        im_spl = np.zeros([num_recs, len(edp_range)])
        for j in range(num_recs):
            for i in range(len(edp_range)):
                if edp_range[i] <= max(edp[j]):
                    im_spl[j][i] = spl[j](edp_range[i])
                else:
                    im_spl[j][i] = im_spl[j][i - 1]

        # Getting the edp-based im quantiles
        im_qtile = np.zeros([len(qtile_range), len(edp_range)])
        for q in range(len(qtile_range)):
            for i in range(len(edp_range)):
                im_qtile[q][i] = np.quantile(im_spl[:, i], qtile_range[q])

        return im_spl, im_qtile

    def splinefit_IDA(self, im, edp, nrecs):

        # Fit the spline and get the values
        spl = {}
        for i in range(nrecs):
            im[i] = np.sort(im[i])
            ind = np.argsort(im[i])
            edp[i] = edp[i][ind]
            ind_mx = int(np.where(edp[i] == max(edp[i]))[0])
            temp_im = im[i]
            temp_edp = edp[i]
            #        if ind_mx < 15:
            #            temp_im = np.delete(im[i],np.arange(ind_mx+1,len(im[i])))
            #            temp_edp = np.delete(edp[i],np.arange(ind_mx+1,len(edp[i])))
            spl[i] = self.cscvn2([temp_edp, temp_im])

        del ind, temp_im, temp_edp

        return spl

    def ida(self, IMpath, dursPath, dtsPath, res_time=10.):
        """
        Postprocess IDA results
        :param IMpath: str                  Path to IM file
        :param dursPath: str                Path to the file containing durations of each record
        :param dtsPath: str                 Path to the file containing time steps of each record
        :param res_time: float              Free vibrations time added to the end of the record
        :return:
        """
        '''
        Current structure of the single pickle file (a very large file I must say )
        Record -> Runs -> EDP type (0=PFA, 1=Displacement, 2=PSD) -> 
        -> array of shape [n x Record length]
        where n is nst for PSD, and nst + 1 for PFA and Displacement
        '''
        """
        The postprocessed file will have the following structure (looks a bit uncomfy though): 
        1. IDA; 2. summary_results (2 keys)
        1.1 ground motions (n_rec keys) -> IM, ISDR, PFA, RISDR (4 keys) -> each being a list of size of number of runs
        2.1 ground motions (n_rec keys) -> IM levels (n_runs keys) -> maxFA, maxISDR (2 keys) -> 
        -> number of storeys (nst for maxISDR) /floors (nst+1 for maxFA) keys -> a single value
        """
        # Read the durations of the records
        durs = list(pd.read_csv(dursPath, header=None)[0])
        # Time steps
        dts = list(pd.read_csv(dtsPath, header=None)[0])
        # Number of records
        nrecs = len(durs)

        # Read the IDA outputs
        if os.path.isdir(self.path):
            data = {}
            for rec in range(nrecs):
                data[rec] = {}

            for filename in os.listdir(self.path):
                if os.path.isfile(self.path / filename):
                    # With no particular order
                    base_filename = os.path.basename(filename).replace(".pickle", "").replace("Record", "").\
                        replace("Run", "").split("_")
                    rec = int(base_filename[0]) - 1
                    run = int(base_filename[1])
                    with open(self.path / filename, 'rb') as file:
                        # idx=0 - accelerations (nst + 1)
                        # idx=1 - displacements (nst + 1)
                        # idx=2 - drifts (nst)
                        data[rec][run] = pickle.load(file)
        else:
            with open(self.path, 'rb') as file:
                data = pickle.load(file)

        # Read the IDA IM levels
        IM = np.genfromtxt(IMpath, delimiter=',')

        # Number of runs per each record
        nruns = len(data[list(data.keys())[0]])

        # Initialize some variables
        im = np.zeros([nrecs, nruns + 1])
        idx = np.zeros([nrecs, nruns], dtype='i')
        mpfa_us = np.full([nrecs, nruns], np.nan)
        mpsd_us = np.full([nrecs, nruns], np.nan)
        mrpsd_us = np.full([nrecs, nruns], np.nan)
        mtdisp_us = np.full([nrecs, nruns + 1], np.nan)
        mtrx = np.full([nrecs, nruns + 1], np.nan)
        mpfa = np.zeros([nrecs, nruns + 1])
        mpsd = np.zeros([nrecs, nruns + 1])
        mtdisp = np.zeros([nrecs, nruns + 1])

        # Initialize target dictionary with its first stage
        res = {'IDA': {}, 'summary_results': {}}
        resKeys = list(res.keys())
        cache = {}
        results = {}

        n_dir = 2 if self.flag3d else 1
        # Loop for each direction for a 3D model
        for d in range(n_dir):
            cache[d] = {}
            results[d] = {}
            # Loop for each record
            for rec in range(1, nrecs+1):
                print("\ngm_%s" % rec)

                # Second stage of the dictionary
                res[resKeys[0]][rec] = {'IM': [], 'ISDR': [], 'PFA': [], 'RISDR': []}
                res[resKeys[1]][rec] = {}

                # Add IM values into the results file
                res[resKeys[0]][rec]["IM"] = IM[rec - 1]
                
                # Sort the IM values
                im[rec - 1, 1:] = np.sort(IM[rec - 1])
                idx[rec - 1, :] = np.argsort(IM[rec - 1])

                # Third stage of the dictionary
                for i in im[rec - 1, 1:]:
                    i = str(np.round(i, 2))
                    res[resKeys[1]][rec][i] = {'maxFA': {}, 'maxISDR': {}, 'maxRISDR': {}}

                # Loop over each run
                for run in range(1, nruns + 1):
                    print(f"gm_{rec}, {run}")
                    # Select analysis results of rec and run
                    selection = data[rec - 1][run]

                    # Get PFAs in g=
                    pfa = np.amax(abs(selection[0][:, 1:]), axis=2)[d]

                    # IML in g
                    iml = str(np.round(IM[rec - 1][run - 1], 2))

                    for st in range(len(pfa)):
                        res[resKeys[1]][rec][iml]["maxFA"][st] = pfa[st]
                    mpfa_us[rec - 1, run - 1] = max(pfa)

                    # Get PSDs in %
                    psd = np.amax(abs(selection[2]), axis=2)[d]

                    for st in range(len(psd)):
                        res[resKeys[1]][rec][iml]["maxISDR"][st + 1] = psd[st]
                    mpsd_us[rec - 1, run - 1] = max(psd)

                    # Getting the residual PSDs in %
                    # Analysis time step
                    dt = dts[rec - 1]
                    idxres = int(durs[rec - 1] / dt)
                    resDrifts = selection[2][:, :, idxres:][d]

                    for st in range(len(psd)):
                        if len(resDrifts[st]) > 0:
                            res[resKeys[1]][rec][iml]["maxRISDR"][st + 1] = sum(resDrifts[st]) / len(resDrifts[st])
                        else:
                            res[resKeys[1]][rec][iml]["maxRISDR"][st + 1] = selection[2][0][st][-1]

                    # Record the peak value of residual drift at each run for each record
                    mrpsd_us[rec - 1, run - 1] = max(np.sum(resDrifts, axis=1) / resDrifts.shape[1])

                    # Get the top displacement in m
                    top_disp = np.amax(abs(selection[1]), axis=n_dir)[d]
                    mtdisp_us[rec - 1, run - 1] = top_disp[-1]

                # Sort the results
                res["IDA"][rec]["PFA"] = mpfa_us[rec - 1, :]
                res["IDA"][rec]["ISDR"] = mpsd_us[rec - 1, :]
                res["IDA"][rec]["RISDR"] = mrpsd_us[rec - 1, :]

                # Repopulate nans with max of data
                # res["IDA"][rec]["RISDR"] = [max(res['IDA'][rec]['RISDR']) if math.isnan(x) else x for
                # x in res['IDA'][rec]['RISDR']]

                mpfa[rec - 1, 1:] = mpfa_us[rec - 1, :][idx[rec - 1]]
                mpsd[rec - 1, 1:] = mpsd_us[rec - 1, :][idx[rec - 1]]
                mtdisp[rec - 1, 1:] = mtdisp_us[rec - 1, :][idx[rec - 1]]

            # Fit the splines to the data
            mtdisp_range = np.linspace(0.01, 1, 200)

            # Quantile ranges to visualize for the IDAs
            qtile_range = np.array([0.16, 0.5, 0.84])

            im_spl = np.zeros([nrecs, len(mtdisp_range)])
            im_spl[:] = np.nan

            # Get the fitted IDA curves for each record
            for rec in range(nrecs):
                interpolator = interp1d(mtdisp[rec], im[rec])

                for i in range(len(mtdisp_range)):
                    if mtdisp_range[i] <= max(mtdisp[rec]):
                        im_spl[rec][i] = interpolator(mtdisp_range[i])
                        if im_spl[rec][i] < im_spl[rec][i - 1]:
                            im_spl[rec][i] = im_spl[rec][i - 1]
                    else:
                        im_spl[rec][i] = im_spl[rec][i - 1]

            # Get the IDA quantiles
            im_qtile = np.zeros([len(qtile_range), len(mtdisp_range)])
            for q in range(len(qtile_range)):
                for i in range(len(mtdisp_range)):
                    im_qtile[q][i] = np.quantile(im_spl[:, i], qtile_range[q])

            # Creating a dictionary for the spline fits
            cache[d] = {"im_spl": im_spl.copy(), "disp": mtdisp.copy(), "im": im.copy(), "im_qtile": im_qtile.copy(),
                        "mtdisp": mtdisp_range.copy()}
            results[d] = res.copy()

        # Exporting
        if self.export:
            if not self.export_path:
                self.export_path = self.path.parents[0]

            self.export_results(self.export_path / "ida_processed", results, "pickle")
            self.export_results(self.export_path / "ida_cache", cache, "pickle")

            print("[SUCCESS] Postprocesssing complete. Results have been exported!")
        else:
            print("[SUCCESS] Postprocesssing complete.")

        return results, cache

    def ida_im_based(self, IMpath, dursPath, res_time=10.):
        """
        Postprocess IDA results
        :param IMpath: str                  Path to IM file
        :param dursPath: str                Path to the file containing durations of each record
        :param res_time: float              Free vibrations time added to the end of the record
        :return:
        """
        '''
        Current structure of the single pickle file (a very large file I must say )
        Record -> Runs -> EDP type (0=PFA, 1=Displacement, 2=PSD) -> 
        -> array of shape [n x Record length]
        where n is nst for PSD, and nst + 1 for PFA and Displacement
        '''
        """
        The postprocessed file will have the following structure (looks a bit uncomfy though): 
        1. IDA; 2. summary_results (2 keys)
        1.1 ground motions (n_rec keys) -> IM, ISDR, PFA, RISDR (4 keys) -> each being a list of size of number of runs
        2.1 ground motions (n_rec keys) -> IM levels (n_runs keys) -> maxFA, maxISDR (2 keys) -> 
        -> number of storeys (nst for maxISDR) /floors (nst+1 for maxFA) keys -> a single value
        """
        # Read the IDA outputs
        with open(self.path, 'rb') as file:
            data = pickle.load(file)

        # Read the IDA IM levels
        IM = np.genfromtxt(IMpath, delimiter=',')

        # Read the durations of the records
        durs = list(pd.read_csv(dursPath, header=None)[0])

        # Number of records
        nrecs = len(data)
        # Number of runs per each record
        nruns = len(data[list(data.keys())[0]])

        # Initialize some variables
        im = np.zeros([nrecs, nruns + 1])
        idx = np.zeros([nrecs, nruns], dtype='i')
        mpfa_us = np.full([nrecs, nruns], np.nan)
        mpsd_us = np.full([nrecs, nruns], np.nan)
        mrpsd_us = np.full([nrecs, nruns], np.nan)
        mtdisp_us = np.full([nrecs, nruns + 1], np.nan)
        mtrx = np.full([nrecs, nruns + 1], np.nan)
        mpfa = np.zeros([nrecs, nruns + 1])
        mpsd = np.zeros([nrecs, nruns + 1])
        mtdisp = np.zeros([nrecs, nruns + 1])

        # Initialize target dictionary with its first stage
        res = {'IDA': {}, 'summary_results': {}}
        resKeys = list(res.keys())

        # Loop for each record
        for rec in range(1, nrecs + 1):
            print("gm_%s" % rec)

            # Second stage of the dictionary
            res[resKeys[0]][rec] = {'IM': [], 'ISDR': [], 'PFA': [], 'RISDR': []}
            res[resKeys[1]][rec] = {}

            # Add IM values into the results file
            res[resKeys[0]][rec]["IM"] = IM[rec - 1]

            # Sort the IM values
            im[rec - 1, 1:] = np.sort(IM[rec - 1])
            idx[rec - 1, :] = np.argsort(IM[rec - 1])

            # Third stage of the dictionary
            for i in im[rec - 1, 1:]:
                i = str(np.round(i, 2))
                res[resKeys[1]][rec][i] = {'maxFA': {}, 'maxISDR': {}, 'maxRISDR': {}}

            # Loop over each run
            for run in range(1, nruns + 1):
                # Select analysis results of rec and run
                selection = data[rec - 1][run]

                # Get PFAs in g
                pfa = np.amax(abs(selection[0][:, 1:]), axis=1)

                # IML in g
                iml = str(np.round(IM[rec - 1][run - 1], 2))

                for st in range(len(pfa)):
                    res[resKeys[1]][rec][iml]["maxFA"][st] = pfa[st]
                mpfa_us[rec - 1, run - 1] = max(pfa)

                # Get PSDs in %
                psd = np.amax(abs(selection[2]), axis=1)

                for st in range(len(psd)):
                    res[resKeys[1]][rec][iml]["maxISDR"][st + 1] = psd[st]
                mpsd_us[rec - 1, run - 1] = max(psd)

                # Getting the residual PSDs
                # Analysis time step
                dt = (durs[rec - 1] + res_time) / selection[0].shape[1]
                idxres = int((durs[rec - 1] - res_time) / dt)

                resDrifts = selection[2][:, idxres:]
                for st in range(len(psd)):
                    res[resKeys[1]][rec][iml]["maxRISDR"][st + 1] = sum(resDrifts[st]) / len(resDrifts[st])
                # Record the peak value of residual drift at each run for each record
                mrpsd_us[rec - 1, run - 1] = max(np.sum(resDrifts, axis=1) / resDrifts.shape[1])

                # Get the top displacement in m
                top_disp = np.amax(abs(selection[1]), axis=1)
                mtdisp_us[rec - 1, run - 1] = top_disp[-1]

            # Sort the results
            res["IDA"][rec]["PFA"] = mpfa_us[rec - 1, :]
            res["IDA"][rec]["ISDR"] = mpsd_us[rec - 1, :]
            res["IDA"][rec]["RISDR"] = mrpsd_us[rec - 1, :]

            # Repopulate nans with max of data
            # res["IDA"][rec]["RISDR"] = [max(res['IDA'][rec]['RISDR']) if math.isnan(x) else x for
            # x in res['IDA'][rec]['RISDR']]

            mpfa[rec - 1, 1:] = mpfa_us[rec - 1, :][idx[rec - 1]]
            mpsd[rec - 1, 1:] = mpsd_us[rec - 1, :][idx[rec - 1]]
            mtdisp[rec - 1, 1:] = mtdisp_us[rec - 1, :][idx[rec - 1]]

        # Fit the splines to the data
        mtdisp_range = np.linspace(0.01, 1, 200)

        # Quantile ranges to visualize for the IDAs
        qtile_range = np.array([0.16, 0.5, 0.84])

        im_spl = np.zeros([nrecs, len(mtdisp_range)])
        im_spl[:] = np.nan

        # Get the fitted IDA curves for each record
        for rec in range(nrecs):
            interpolator = interp1d(mtdisp[rec], im[rec])

            for i in range(len(mtdisp_range)):
                if mtdisp_range[i] <= max(mtdisp[rec]):
                    im_spl[rec][i] = interpolator(mtdisp_range[i])
                    if im_spl[rec][i] < im_spl[rec][i - 1]:
                        im_spl[rec][i] = im_spl[rec][i - 1]
                else:
                    im_spl[rec][i] = im_spl[rec][i - 1]

        # Get the IDA quantiles
        im_qtile = np.zeros([len(qtile_range), len(mtdisp_range)])
        for q in range(len(qtile_range)):
            for i in range(len(mtdisp_range)):
                im_qtile[q][i] = np.quantile(im_spl[:, i], qtile_range[q])

        # Similarly get the fitted IDA curves for each record on IM-basis
        # Sort the displacements
        disp_sorted = np.sort(mtdisp)
        idx_disp = np.argsort(mtdisp)
        # Get the corresponding IM values
        im_disp = np.array(list(map(lambda x, y: y[x], idx_disp, im)))
        # IML range to use for fitting
        im_range = np.linspace(0.0, np.max(im), 50)

        disp_spl = np.zeros([nrecs, len(im_range)])
        disp_spl[:] = np.nan
        for rec in range(nrecs):
            interpolator = interp1d(im_disp[rec], disp_sorted[rec], fill_value=max(disp_sorted[rec]),
                                    bounds_error=False)
            disp_spl[rec] = interpolator(im_range)

        # Creating a dictionary for the spline fits
        cache = {"im_spl": im_spl, "disp": mtdisp, "im": im, "im_qtile": im_qtile, "mtdisp": mtdisp_range}

        # Exporting
        if self.export:
            self.export_results(self.path.parents[0] / "ida_processed", res, "pickle")
            self.export_results(self.path.parents[0] / "ida_cache", cache, "pickle")
            print("[SUCCESS] Postprocesssing complete. Results have been exported!")
        else:
            print("[SUCCESS] Postprocesssing complete.")

        return res, cache

    def mafe_direct_im_based(self, eta, beta, sa_haz, Hs):
        """
        Details:
        Compute the MAFE of a limit state defined via a fitted lognormal
        distribution by integrating directly with the  hazard curve
        Treat the hazard input to avoid errors.
        We strip out:
         1. the negative H values (usually at the beginning)
         2. the points with constant s (usually at the end)
        Information:
        Author: Gerard J. O'Reilly
        First Version: April 2020
        Notes:
        References:
        Porter KA, Beck JL, Shaikhutdinov R V. Simplified Estimation of Economic
        Seismic Risk for Buildings. Earthquake Spectra 2004; 20(4):
        1239–1263. DOI: 10.1193/1.1809129.
        Inputs:
        :param eta: float                           Fragility function median (intensity)
        :param beta: float                          Fragility function dispersion (total)
        :param sa_haz: array                        List of corresponding intensities for Hs
        :param Hs: array                            List of values of annual probability of exceedance
        :return: float                              Mean annual frequency of exceedance
        """

        # Do first strip
        s_f = []
        H_f = []
        for aa, bb in zip(sa_haz, Hs):
            if bb > 0:
                s_f.append(aa)
                H_f.append(bb)

        # Do second strip
        s_ff = []
        H_ff = []
        for i in range(len(s_f) - 1):
            if H_f[i] - Hs[i + 1] > 0:
                s_ff.append(s_f[i])
                H_ff.append(H_f[i])
        s_ff.append(s_f[-1])
        H_ff.append(H_f[-1])

        # Overwrite the initial variable for convenience
        s = s_ff
        H = H_ff

        # First we compute the PDF value of the fragility at each of the discrete
        # hazard curve points
        p = lognorm.cdf(s, beta, scale=eta)

        # This function computes the MAF using Method 1 outlined in
        # Porter et al. [2004]
        # This assumes that the hazard curve is linear in logspace between
        # discrete points among others

        # Initialise some arrays
        ds = []
        ms = []
        dHds = []
        dp = []
        dl = []

        for i in np.arange(len(s) - 1):
            ds.append(s[i + 1] - s[i])
            ms.append(s[i] + ds[i] * 0.5)
            dHds.append(np.log(H[i + 1] / H[i]) / ds[i])
            dp.append(p[i + 1] - p[i])
            dl.append(p[i] * H[i] * (1 - np.exp(dHds[i] * ds[i])) - dp[i] / ds[i] * H[i] * (
                        np.exp(dHds[i] * ds[i]) * (ds[i] - 1 / dHds[i]) + 1 / dHds[i]))

        # Compute the MAFE
        l = sum(dl)
        return l

    def verify_mafc(self, res, hazardPath, targetMAFC, MApath, ipbsdPath, period=None):
        """
        Verifies that MAFC is below the target value
        :param res: pickle
        :param hazardPath: str
        :param targetMAFC: float
        :param MApath: str
        :param ipbsdPath: str
        :param period: float
        :return: bool
        """
        # Read the hazard information
        with open(hazardPath, 'rb') as file:
            [im, s, apoe] = pickle.load(file)
            im = np.array(im)
            s = np.array(s)
            apoe = np.array(apoe)

        # Read the fundamental period of the structure
        if period is None:
            with open(MApath) as f:
                results = json.load(f)
                period = results["Periods"][0]

        # IPBSD results
        with open(ipbsdPath, "rb") as f:
            ipbsd = pickle.load(f)
            gamma = ipbsd["part_factor"]

        # IDA results from the NTHA
        im_spl = res["im_spl"]
        disp_range = res["mtdisp"]
        # Interpolation functions
        spl_interp = interp1d(disp_range, im_spl*gamma)

        # Top displacement capacity at collapse (used to find the collapse IM distribution)
        cap_disp = disp_range[-1]
        # Collapse IM distribution and sorted
        spl_mu = spl_interp(cap_disp)
        spl_mu = np.sort(spl_mu)

        # Median and standard deviation of the collapse IM distribution
        eta = np.median(spl_mu)
        beta = np.std(np.log(spl_mu))
        print(f"Collapse capacity: eta={eta:.2f}, beta_RTR={beta:.2f}, period={period:.2f}")

        indx_T = int(period * 10)

        l = self.mafe_direct_im_based(eta, beta, s[indx_T], apoe[indx_T])
        if l <= targetMAFC:
            print(f"[SUCCESS] Actual MAFC over target MAFC: {l/targetMAFC:.2f}. Actual MAFC: {l:.5f}")
            return True
        else:
            print(f"[FAILURE] Actual MAFC over target MAFC: {l/targetMAFC:.2f}. Actual MAFC: {l:.5f}")
            return False

