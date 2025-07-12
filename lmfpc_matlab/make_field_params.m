function field_params = make_field_params(field_choice, particle_params)

    mime = particle_params.mime;
    tite = particle_params.tite;

    % Default beta
    bi = NaN; be = NaN; waveT = NaN; kvalue = NaN;

    switch field_choice
        case {-41, -45, -47, -451}
            bi = 1.; waveT = 110.72; kvalue = 0.05;
        case {-49, -491}
            bi = 1.; waveT = 110.72; kvalue = -0.05;
        case {-495, -4951}
            bi = 1.; waveT = 110.72; kvalue = [-0.05, 0.05];
        case -43
            bi = 1.; waveT = 55.61; kvalue = 0.1;
        case {-51, -53, -511}
            bi = 1.; waveT = 20.1657; kvalue = 0.8;
        case -55
            bi = 1.; waveT = 25.78599; kvalue = 0.4;
        case -59
            bi = 1.; waveT = 22.44201; kvalue = 0.6;
        case -61
            bi = 1.; waveT = 23.94790; kvalue = 0.5;
        case {-63, -6300}
            bi = 1.; waveT = 23.57534; kvalue = 0.525;
        case -65
            bi = 1.; waveT = 23.57534; kvalue = -0.525;
        case {-653, -6531, -6532, -6533}
            bi = 1.; waveT = 23.57534; kvalue = [-0.525, 0.525];
        case -6301
            bi = 0.1; waveT = 10.61342774; kvalue = 0.525;
        case -6303
            bi = 0.3; waveT = 14.52614851; kvalue = 0.525;
        case -633
            bi = 3.; waveT = 41.67231726; kvalue = 0.525;
        case -6310
            bi = 10.; waveT = 90.66208508; kvalue = 0.525;
        otherwise
            error('Unsupported field_choice: %d', field_choice);
    end

    be = bi * mime / tite;
    em_eps = 0.02;  
    delta_phi = 0.;

    field_params = struct( ...
        'field_choice', field_choice, ...
        'bi', bi, ...
        'be', be, ...
        'waveT', waveT, ...
        'kvalue', kvalue, ...
        'em_eps', em_eps, ...
        'delta_phi', delta_phi ...
    );


