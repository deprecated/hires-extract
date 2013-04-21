import image_utils

p85fac = 500.0/600.0            # exposure time ratio p85/p84


def line_stages_two_slit(id0, id1, drange, pmax, psmooth):
    image_utils.line_stages(
        "p84", id0, id1, drange=drange, pmax=pmax, psmooth=psmooth)
    image_utils.line_stages(
        "p85", id0, id1, drange=[x*p85fac for x in drange],
        pmax=pmax*p85fac, psmooth=2.0)


line_stages_two_slit("Fe_II_5159", "[Fe II] 5159",
                     drange=[-2.0, 60.0], pmax=40.0, psmooth=2.0)
line_stages_two_slit("He_I_T_5876", "He I 5876",
                     drange=[-5.0, 1000.0], pmax=600.0)
line_stages_two_slit("N_II_5755", "[N II] 5755",
                     drange=[-10.0, 800.0], pmax=600.0)
line_stages_two_slit("N_II_6583", "[N II] 6583",
                     drange=[-10.0, 60000.0], pmax=32000.0)
line_stages_two_slit("N_II_6548", "[N II] 6548",
                     drange=[-10.0, 15000.0], pmax=8000.0)
line_stages_two_slit("N_I_5200", "[N I] 5200",
                     drange=[-2.0, 100.0], pmax=100.0, sky=1)
line_stages_two_slit("O_I_6046", "O I 6046", doublet=True,
                     drange=[-3.0, 180.0], pmax=400.0)
line_stages_two_slit("O_I_7002", "O I 7002", doublet=True,
                     drange=[-3.0, 180.0], pmax=400.0)
line_stages_two_slit("Si_II_6347", "Si II 6347",
                     drange=[-3.0, 160.0], pmax=120.0, psmooth=1.0)
line_stages_two_slit("Si_II_6371", "Si II 6371",
                     drange=[-3.0, 120.0], pmax=80.0, psmooth=1.0)
line_stages_two_slit("N_I_5198", "[N I] 5198",
                     drange=[-2.0, 140.0], pmax=160.0, sky=1)
line_stages_two_slit("O_I_6300", "[O I] 6300",
                     drange=[-10.0, 1600.0], sky=1, pmax=1600.0)
line_stages_two_slit("C_II_6578", "C II 6578",
                     drange=[-2.0, 150.0], pmax=80.0, psmooth=2.0)
line_stages_two_slit("He_I_S_6678", "He I 6678",
                     drange=[-5.0, 700.0], pmax=400.0)
line_stages_two_slit("O_I_5577", "[O I] 5577",
                     drange=[-3.0, 250.0], sky=1, pmax=160.0)
line_stages_two_slit("O_III_5007", "[O III] 5007",
                     drange=[-10.0, 10000.0], pmax=4000.0)
line_stages_two_slit("S_II_6716", "[S II] 6716",
                     drange=[-10.0, 3000.0], pmax=1600.0)
line_stages_two_slit("S_II_6731", "[S II] 6731",
                     drange=[-10.0, 5000.0], pmax=3200.0)
line_stages_two_slit("S_III_6312", "[S III] 6312",
                     drange=[-4.0, 330.0], pmax=320.0)
line_stages_two_slit("Cl_III_5518", "[Cl III] 5518",
                     drange=[-2.0, 50.0], pmax=40.0, psmooth=2.0)
line_stages_two_slit("Cl_III_5538", "[Cl III] 5538",
                     drange=[-2.0, 100.0], pmax=40.0, psmooth=2.0)
line_stages_two_slit("Fe_III_5270", "[Fe III] 5270",
                     drange=[-2.0, 40.0], pmax=20.0, psmooth=2.0)
line_stages_two_slit("Fe_III_4881", "[Fe III] 4881",
                     drange=[-2.0, 60.0], pmax=40.0, psmooth=2.0)
line_stages_two_slit("Fe_II_5262", "[Fe II] 5262",
                     drange=[-2.0, 25.0], pmax=40.0, psmooth=2.0)

