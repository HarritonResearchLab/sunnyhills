def test_env(tls_or_wotan): 

    if tls_or_wotan.lower()=='wotan': 
        import wotan 
        
    else: 

        import transitleastsquares 
        
