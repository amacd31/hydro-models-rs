pub mod hydromodels {
    use std::collections::HashMap;

    /*
    Unit hydrograph ordinates for UH1 derived from S-curves.
    */
    pub(crate) fn s_curves1(t: f64, x4: f64) -> f64 {

        if t <= 0.0 {
            return 0.0;
        } else if t < x4 {
            return (t / x4).powf(2.5);
        } else {
            // t >= x4
            return 1.0;
        }
    }

    /*
    Unit hydrograph ordinates for UH2 derived from S-curves.
    */
    pub(crate) fn s_curves2(t: f64, x4: f64) -> f64 {

        if t <= 0.0 {
            return 0.0;
        } else if t < x4 {
            return 0.5 * (t / x4).powf(2.5);
        } else if t < 2.0 * x4 {
            return 1.0 - 0.5 * (2.0 - t / x4).powf(2.5);
        } else {
            // t >= x4
            return 1.0;
        }
    }

    /*
    Generated simulated streamflow for given rainfall and potential evaporation.

    :param precip: Catchment average rainfall.
    :param potential_evap: Catchment average potential evapotranspiration.
    :param params: x parameters for the model.

    :return: Vector of simulated streamflow.
    */
    pub fn gr4j(precip: &[f64], potential_evap: &[f64], params: HashMap<&str, f64>, states: Option<HashMap<&str, f64>>) -> Vec<f64> {
        let mut qsim: Vec<f64> = Vec::new();

        let x1 = params["X1"];
        let x2 = params["X2"];
        let x3 = params["X3"];
        let x4 = params["X4"];

        let n_uh1 = x4.ceil() as i32;
        let n_uh2 = (2.0 * x4).ceil() as i32;

        let mut uh1_ordinates: Vec<f64> = Vec::with_capacity(n_uh1 as usize);
        for _i in 0..n_uh1 {
            uh1_ordinates.push(0.);
        }

        let mut uh2_ordinates: Vec<f64> = Vec::with_capacity(n_uh2 as usize);
        for _i in 0..n_uh2 {
            uh2_ordinates.push(0.);
        }

        //UH1 = states.get('UH1', [0] * n_uh1)
        let mut uh1: Vec<f64> = Vec::with_capacity(n_uh1 as usize);
        for _i in 0..n_uh1 {
            uh1.push(0.);
        }

        //UH2 = states.get('UH2', [0] * n_uh2)
        let mut uh2: Vec<f64> = Vec::with_capacity(n_uh2 as usize);
        for _i in 0..n_uh2 {
            uh2.push(0.);
        }

        for t in 1..(n_uh1 + 1) {
            uh1_ordinates[(t - 1) as usize] = s_curves1(t as f64, x4) -
                s_curves1(t as f64 - 1., x4);
        }

        for t in 1..(n_uh2 + 1) {
            uh2_ordinates[(t - 1) as usize] = s_curves2(t as f64, x4) -
                s_curves2(t as f64 - 1., x4);
        }

        let mut production_store; // S
        let mut routing_store; // R
        if states.is_some() {
            let states_hash = states.unwrap();
            production_store = *states_hash.get("production_store").unwrap_or(&0.); // S
            routing_store = *states_hash.get("routing_store").unwrap_or(&0.); // R

        } else {
            production_store = 0.; // S
            routing_store = 0.; // R
        }

        for (p, e) in precip.iter().zip(potential_evap) {
            let net_evap;
            let mut routing_pattern;
            let reservoir_production;
            if p > e {
                net_evap = 0.;
                let mut scaled_net_precip = (p - e) / x1;
                if scaled_net_precip > 13. {
                    scaled_net_precip = 13.;
                }
                let tanh_scaled_net_precip = scaled_net_precip.tanh();
                reservoir_production = (x1 * (1. - (production_store / x1).powf(2.)) *
                                            tanh_scaled_net_precip) /
                    (1. + production_store / x1 * tanh_scaled_net_precip);

                routing_pattern = p - e - reservoir_production;
            } else {
                let mut scaled_net_evap = (e - p) / x1;
                if scaled_net_evap > 13. {
                    scaled_net_evap = 13.;
                }
                let tanh_scaled_net_evap = scaled_net_evap.tanh();

                let ps_div_x1 = (2. - production_store / x1) * tanh_scaled_net_evap;
                net_evap = production_store * (ps_div_x1) /
                    (1. + (1. - production_store / x1) * tanh_scaled_net_evap);

                reservoir_production = 0.;
                routing_pattern = 0.;
            }

            production_store = production_store - net_evap + reservoir_production;

            let percolation = production_store /
                (1. + (production_store / 2.25 / x1).powf(4.)).powf(0.25);
            routing_pattern = routing_pattern + (production_store - percolation);

            production_store = percolation;

            for i in 0..(uh1.len() - 1) {
                uh1[i] = uh1[i + 1] + uh1_ordinates[i] * routing_pattern;
            }
            if let (Some(last_uh1), Some(last_ordinate)) = (uh1.last_mut(), uh1_ordinates.last()) {
                *last_uh1 = *last_ordinate * routing_pattern;
            }

            for j in 0..(uh2.len() - 1) {
                uh2[j] = uh2[j + 1] + uh2_ordinates[j] * routing_pattern
            }
            if let (Some(last_uh2), Some(last_ordinate)) = (uh2.last_mut(), uh2_ordinates.last()) {
                *last_uh2 = *last_ordinate * routing_pattern;
            }

            let groundwater_exchange = x2 * (routing_store / x3).powf(3.5);
            routing_store = (0.0f64).max(routing_store + uh1[0] * 0.9 + groundwater_exchange);

            let r2 = routing_store / (1. + (routing_store / x3).powf(4.)).powf(0.25);
            let qr = routing_store - r2;
            routing_store = r2;
            let qd = (0.0f64).max(uh2[0] * 0.1 + groundwater_exchange);
            let q = qr + qd;

            qsim.push(q);
        }

        qsim
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use hydromodels;

    #[test]
    fn gr4j_simple_test() {
        let mut params = HashMap::new();

        params.insert("X1", 10.);
        params.insert("X2", 5.);
        params.insert("X3", 4.);
        params.insert("X4", 1.);

        let mut expected: Vec<f64> = Vec::new();
        expected.push(0.13030636843636356);
        expected.push(0.7348755690383941);
        expected.push(3.8482364547176102);
        expected.push(11.246676991863168);
        expected.push(13.90322269162079);

        let result = hydromodels::gr4j(&[10., 2., 3., 4., 5.], &[0.5, 0.5, 0.5, 0.5, 0.5], params, None);

        assert_eq!(result, expected);
    }

    #[test]
    fn gr4j_scaled_test() {
        let mut params = HashMap::new();

        params.insert("X1", 10.);
        params.insert("X2", 5.);
        params.insert("X3", 4.);
        params.insert("X4", 1.);

        let mut expected: Vec<f64> = Vec::new();

        expected.push(0.13030636843636356);
        expected.push(0.5196067786804051);
        expected.push(0.7978026532390882);
        expected.push(128.2406932741725);
        expected.push(20.992238476264593);

        let result =
            hydromodels::gr4j(&[10., 2., 3., 150., 5.], &[0.5, 14., 0.5, 10., 0.5], params, None);

        assert_eq!(result, expected);
    }

    #[test]
    fn gr4j_set_production_store() {
        let mut params = HashMap::new();

        params.insert("X1", 10.);
        params.insert("X2", 5.);
        params.insert("X3", 4.);
        params.insert("X4", 1.);

        let mut states = HashMap::new();

        states.insert("production_store", 10.);

        let mut expected: Vec<f64> = Vec::new();

        expected.push(5.1602235324393675);
        expected.push(10.091188499285725);
        expected.push(9.82974398339987);
        expected.push(136.95699908376093);
        expected.push(21.019904684254975);

        let result =
            hydromodels::gr4j(&[10., 2., 3., 150., 5.], &[0.5, 14., 0.5, 10., 0.5], params, Some(states));

        assert_eq!(result, expected);
    }

    #[test]
    fn s_curves1_test() {
        assert_eq!(hydromodels::s_curves1(-1.0, 2.0), 0.0);
        assert_eq!(hydromodels::s_curves1(0.0, 2.0), 0.0);
        assert_eq!(hydromodels::s_curves1(1.0, 2.0), 0.1767766952966369);
        assert_eq!(hydromodels::s_curves1(2.0, 2.0), 1.0);
        assert_eq!(hydromodels::s_curves1(2.1, 2.0), 1.0);
    }

    #[test]
    fn s_curves2_test() {
        assert_eq!(hydromodels::s_curves2(-1.0, 2.0), 0.0);
        assert_eq!(hydromodels::s_curves2(0.0, 2.0), 0.0);
        assert_eq!(hydromodels::s_curves2(1.0, 2.0), 0.08838834764831845);
        assert_eq!(hydromodels::s_curves2(2.0, 2.0), 0.5);
        assert_eq!(hydromodels::s_curves2(2.1, 2.0), 0.5601759051904955);
        assert_eq!(hydromodels::s_curves2(4.0, 2.0), 1.0);
        assert_eq!(hydromodels::s_curves2(4.1, 2.0), 1.0);
    }
}
