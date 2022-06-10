def tls_intransit_stats(t, y, transit_times, transit_duration_in_days):
    """
    Return all intransit odd and even flux points
    
    Source: Hippke's TLS Code
    Status: light modification 
    """

    import numpy as np

    all_flux_intransit_odd = np.array([])
    all_flux_intransit_even = np.array([])
    all_time_intransit_odd = np.array([])
    all_time_intransit_even = np.array([])
    all_flux_intransit = np.array([])
    all_idx_intransit = np.array([])
    per_transit_count = np.zeros([len(transit_times)])
    transit_depths = np.zeros([len(transit_times)])
    transit_depths_uncertainties = np.zeros([len(transit_times)])

    for i in range(len(transit_times)):

        depth_mean_odd = np.nan
        depth_mean_even = np.nan
        depth_mean_odd_std = np.nan
        depth_mean_even_std = np.nan

        mid_transit = transit_times[i]
        tmin = mid_transit - 0.5 * transit_duration_in_days
        tmax = mid_transit + 0.5 * transit_duration_in_days
        if np.isnan(tmin) or np.isnan(tmax):
            idx_intransit = []
            flux_intransit = []
            time_intransit = []
            mean_flux = np.nan
        else:
            idx_intransit = np.where(np.logical_and(t > tmin, t < tmax))[0]
            flux_intransit = y[idx_intransit]
            time_intransit = t[idx_intransit]
            if len(y[idx_intransit]) > 0:
                mean_flux = np.mean(y[idx_intransit])
            else:
                mean_flux = np.nan
        intransit_points = np.size(y[idx_intransit])
        transit_depths[i] = mean_flux
        if len(y[idx_intransit] > 0):
            transit_depths_uncertainties[i] = np.std(y[idx_intransit]) / np.sqrt(
                intransit_points
            )
        else:
            transit_depths_uncertainties[i] = np.nan
        per_transit_count[i] = intransit_points

        # Check if transit odd/even to collect the flux for the mean calculations
        if i % 2 == 0:  # even
            all_flux_intransit_even = np.append(
                all_flux_intransit_even, flux_intransit
            )

            all_time_intransit_even = np.append(
                all_time_intransit_even, time_intransit
            )
        else:  # odd
            all_flux_intransit_odd = np.append(
                all_flux_intransit_odd, flux_intransit
            )
            all_time_intransit_odd = np.append(
                all_time_intransit_odd, time_intransit
            )

        if len(all_flux_intransit_odd) > 0:
            depth_mean_odd = np.mean(all_flux_intransit_odd)

            depth_mean_odd_std = np.std(all_flux_intransit_odd) / np.sum(
                len(all_flux_intransit_odd)
            ) ** (0.5)
        if len(all_flux_intransit_even) > 0:
            depth_mean_even = np.mean(all_flux_intransit_even)
            depth_mean_even_std = np.std(all_flux_intransit_even) / np.sum(
                len(all_flux_intransit_even)
            ) ** (0.5)

    return (
        depth_mean_odd,
        depth_mean_even,
        depth_mean_odd_std,
        depth_mean_even_std,
        (all_time_intransit_odd,all_flux_intransit_odd),
        (all_time_intransit_even, all_flux_intransit_even),
        per_transit_count,
        transit_depths,
        transit_depths_uncertainties,
    )

