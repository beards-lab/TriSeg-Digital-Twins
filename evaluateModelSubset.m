function cost = evaluateModelSubset(m_opt, m_initial, idx, UWpatients, PatID)
    m_full = m_initial;
    m_full(idx) = m_opt;

    cost = evaluateModelUW(m_full, UWpatients, PatID);
end