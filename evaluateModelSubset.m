function cost = evaluateModelSubset(m_opt, m_initial, idx, UWpatients, PATIENT_NO,Geo_Opt)
    m_full = m_initial;
    m_full(idx) = m_opt;
    cost = evaluateModelUW(m_full, UWpatients, PATIENT_NO,Geo_Opt);
end