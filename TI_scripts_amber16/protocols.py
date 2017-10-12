"""
Template of the amber 16 inputs
"""
STD_INPUT = {'min': """minimisation
                         &cntrl
                           imin = 1, ntmin = 2, maxcyc = 8000,
                           ntpr = 20, ntwe = 20,
                           dx0 = 1.0D-7,
                           ntb = 1, barostat = 2,

                           icfe = 1, ifsc = 1, clambda = {clambda:.3f}, scalpha = 0.5, scbeta = 12.0,
                           logdvdl = 0,
                           {timask1} {timask2}
                           {scmask1} {scmask2}

                         /

                         &ewald
                         / """,

            'heat':"""heating
                     &cntrl
                       imin = 0, nstlim = 80000, irest = 0, ntx = 1, dt = 0.002,
                       ntt = 3, temp0 = 300.0, tempi = 50.0, tautp = 1.0,
                       ntc = 2, ntf = 1,  barostat = 2,
                       ntb = 1,
                       ioutfm = 1, iwrap = 1,
                       ntwe = 1000, ntwx = 1000, ntpr = 1000, ntwr = 5000,

                       nmropt = 1,
                       ntr = 1, restraint_wt = 5.00,
                       restraintmask='!:WAT | !@H=',

                       icfe = 1, ifsc = 1, clambda = {clambda:.3f}, scalpha = 0.5, scbeta = 12.0,
                       logdvdl = 0,
                       {timask1} {timask2}
                       {scmask1} {scmask2}


                     /

                     &ewald
                     /

                     &wt
                       type='TEMP0',
                       istep1 = 0, istep2 = 8000,
                       value1 = 50.0, value2 = 300.0
                     /

                     &wt type = 'END'
                     /
                    """,
            'prod':"""TI simulation
                     &cntrl
                       imin = 0, nstlim = {psteps}, irest = 1, ntx = 5, dt = 0.002,
                       ntt = 3, temp0 = 300.0, gamma_ln = 2.0, ig = -1,
                       ntc = 2, ntf = 1,
                       ntb = 2, barostat = 2,
                       ntp = 1, pres0 = 1.0, taup = 2.0,
                       ioutfm = 1, iwrap = 1,
                       ntwe = 1000, ntwx = 10000, ntpr = 10000, ntwr = 20000,

                       icfe = 1, ifsc = 1, clambda = {clambda:.3f}, scalpha = 0.5, scbeta = 12.0,
                       logdvdl = 1,
                       ifmbar = 0, bar_intervall = 1000, bar_l_min = 0.0, bar_l_max = 1.0,
                         bar_l_incr = 0.1,
                            {timask1} {timask2}
                            {scmask1} {scmask2}


                     /


                     &ewald
                     / """}


GEN_INPUT = {'1-min':"""&cntrl
                            imin=1,
                            maxcyc=5000,
                            ncyc=200,
                            drms=0.5,
                            ntr=1,
                            ntpr = 10, ntwe = 10,
                            ntmin = 2,
                            dx0 = 1.0D-7,
                            ntb = 1,
                            barostat = 2,
                            restraintmask="!:WAT,NA,CL,Na+,Cl-,IP,IM",
                            restraint_wt=50,
                            icfe = 1, ifsc = 1, clambda = {clambda:.3f}, scalpha = 0.5, scbeta = 12.0,
                            logdvdl = 0,
                            {timask1} {timask2}
                            {scmask1} {scmask2}
                            /
                            &ewald
                            /
                            &end
                            """,

              '2-min':"""minimisation
                        &cntrl
                        imin = 1, ntmin = 2, maxcyc = 5000, ncyc=500,
                        ntpr = 1000, ntwe = 1000,
                        dx0 = 1.0D-7,
                        ntb = 1,
                        barostat = 2,
                        restraintmask="@CA,C,N,O",
                        restraint_wt=1,
                        icfe = 1, ifsc = 1, clambda = {clambda:.3f}, scalpha = 0.5, scbeta = 12.0,
                        logdvdl = 0,
                        {timask1} {timask2}
                        {scmask1} {scmask2}

                        /
                        &ewald
                        /
                        &end""",
              '3-equil': """heating
                            &cntrl
                            imin = 0, nstlim = 200000, irest = 0, ntx = 1, dt = 0.001,
                            ntt = 3, tautp = 1.0, temp0=300.0,tempi=50.0,
                            ntc = 2, ntf = 1,
                            ntb = 1,
                            barostat = 2,
                            ioutfm = 1, iwrap = 0,
                            ntwe = 1000, ntwx = 1000, ntpr = 1000, ntwr = 1000,

                            nmropt = 1,
                            ntr = 1, restraint_wt = 10.00,
                            restraintmask='!:WAT',

                            icfe = 1, ifsc = 1, clambda = {clambda:.3f}, scalpha = 0.5, scbeta = 12.0,
                            logdvdl = 0,
                            {timask1} {timask2}
                            {scmask1} {scmask2}
                            /
                            &ewald
                            /
                            &wt
                            type='TEMP0',
                            istep1 = 0, istep2 = 100000,
                            value1 = 50.0, value2 = 300.0
                            /
                            &wt type = 'END'
                            /""",

              '4-equil': """Equilibrarion
                        &cntrl
                        imin = 0, nstlim = 200000, irest = 1, ntx = 5, dt = 0.002,
                        ntt = 3, temp0 = 300.0, gamma_ln = 1.0, ig = -1,
                        ntc = 2, ntf = 1,
                        ntb = 2,
                        ntp = 1, pres0 = 1.0, taup = 2.0,
                        ioutfm = 1, iwrap = 0,
                        ntwe = 1000, ntwx = 1000, ntpr = 1000, ntwr = 1000,

                        icfe = 1, ifsc = 1, clambda = {clambda:.3f}, scalpha = 0.5, scbeta = 12.0,
                        logdvdl = 1,
                        barostat = 2, ifmbar = 0, bar_intervall = 1000, bar_l_min = 0.0, bar_l_max = 1.0,
                        bar_l_incr = {increment},
                        {timask1} {timask2}
                        {scmask1} {scmask2}

                        /


                        &ewald
                        / """,


              'prod' : """TI simulation
                        &cntrl
                        imin = 0, nstlim = {psteps}, irest = 1, ntx = 5, dt = 0.002,
                        ntt = 3, temp0 = 300.0, gamma_ln = 2.0, ig = -1,
                        ntc = 2, ntf = 1,
                        ntb = 2,
                        ntp = 1, pres0 = 1.0, taup = 2.0,
                        ioutfm = 1, iwrap = 0,
                        ntwe = 1000, ntwx = 10000, ntpr = 10000, ntwr = 10000,

                        icfe = 1, ifsc = 1, clambda = {clambda:.3f}, scalpha = 0.5, scbeta = 12.0,
                        logdvdl = 1,
                        barostat = 2, ifmbar = 0, bar_intervall = 1000, bar_l_min = 0.0, bar_l_max = 1.0,
                        bar_l_incr = {increment},
                        {timask1} {timask2}
                        {scmask1} {scmask2}

                        /


                        &ewald
                        / """}


            #   '4-equil': """Equilibrarion
            #                 &cntrl
            #                 imin = 0, nstlim = 250000, irest = 0, ntx = 1, dt = 0.002,
            #                 ntt = 3, temp0 = 300.0, tempi = 150.0, tautp = 1.0,
            #                 ntc = 2, ntf = 2,
            #                 ntb = 1, pres0 = 1.0,
            #                 ioutfm = 1, iwrap = 1,
            #                 ntwe = 1000, ntwx = 1000, ntpr = 1000, ntwr = 5000,
              #
            #                 nmropt = 1,
            #                 ntr = 1, restraint_wt = 10.00,
            #                 restraintmask='@CA,C,N,O',
              #
            #                 icfe = 1, ifsc = 1, clambda = {clambda:.3f}, scalpha = 0.5, scbeta = 12.0,
            #                 logdvdl = 0,
            #                 {timask1} {timask2}
            #                 {scmask1} {scmask2}
            #                 /
            #                 &ewald
            #                 &wt type = 'END'
            #                 /""",
