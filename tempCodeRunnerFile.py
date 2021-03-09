e(n): 
        fx_i = p[i][1]
        
        # forward diff
        if (i + h_dir) == n: 
            print(f'[F] f\'({p[i][0]}) = can\'t be calculated using forward m since f(x_i + h) wasn\'t provided.')
        else:    
            fx_i_ph = p[i + h_dir][1]
            df_f = ( fx_i_ph - fx_i ) / h
            p[i][2] = df_f            
            print(f'[F] f\'({p[i][0]}) = {df_f:<1.5f}')

        # backward diff
        if (i - h_dir) == -1:
            print(f'[B] f\'({p[i][0]}) = can\'t be calculated using backward m since f(x_i - h) wasn\'t provided.')
        else:
            fx_i_mh = p[i - h_dir][1]
            df_b = ( fx_i - fx_i_mh ) / h
            p[i][3] = df_b
            print(f'[B] 