import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
import numpy as np



def visualize_sub(pd_net, targets, name):
    tf_genes = ['Mef2a', 'Mef2c', 'Nr2f2', 'Tbx5', 'Gata4', 'Meis2', 'Plagl1', 'Isl1', 'Pitx2', 'Tbx20', 'Hand2', 'Nkx2-5']
    
    targets = targets[:100]

    string = ''
    for s in targets:
        string += f'{s}, '

    targets += tf_genes

    net, G1, G2 = plot_subnetwork(pd_net, tf_genes, targets)
    # net = plot_subnetwork(links_ko_e85.filtered_links['IFT-CMs_KO'], tf_genes, targets)
    # net.toggle_physics(False)

    # net.barnes_hut()
    #net.show('net.html')
    net.show_buttons(filter_=['physics'])
    net.show(name, notebook=False)


def plot_subnetwork(network, tf_genes, targets):
    graph_df = {'source': [],
                 'target': [],
                 'weight': [],
                 'edge_cols': []}

    for curr_gene in tf_genes:
        curr_gene_df = network[network['source'] == curr_gene]
        gene_df_targs = curr_gene_df[curr_gene_df['target'].isin(targets)]

        for index, row in gene_df_targs.iterrows():
            graph_df['source'].append(row['source'])
            graph_df['target'].append(row['target'])
            graph_df['weight'].append(row['coef_mean'])
            if row['coef_mean'] > 0:
                graph_df['edge_cols'].append('green')
            else:
                graph_df['edge_cols'].append('red')

    graph_df = pd.DataFrame(graph_df)

    G = nx.from_pandas_edgelist(graph_df, source='source',
                                    target='target',
                                    edge_attr='weight',
                                    create_using=nx.DiGraph())

    from pyvis.network import Network
    net = Network(notebook=True, directed=True, height='500px', width='1000px')
    net.from_nx(G)

    for edge in net.edges:
        if edge['width'] > 0:
            edge['color'] = 'grey'
            edge['arrowStrikethrough'] = False
        else:
            edge['color'] = 'red'
            edge['arrows'] = {"to": {"enabled": True, "type": "bar",
                                     "scaleFactor":1.4}}

        edge['value'] = np.abs(edge['width'])


    for node in net.nodes:
        if node['id'] in tf_genes:
            #node['shape'] = 'circle'
            node['color'] = 'darkgrey'
            node['size'] = 10
            node['font'] = {'size': 50,
                            'color': 'black',
                            'style': 'italic',
                            'multi': 'html'}
        else:
            node['color'] = 'darkgrey'
            node['font'] = {'size': 26,
                            'color': 'black',
                            'style':'italic',
                            'multi': 'html'}
            node['size'] = 6 

        if node['id'] in ['Gata4', 'Nr2f2']:
            node['shape'] = 'diamond'
            node['color'] = 'darkgrey'
            node['size'] = 20
            node['font'] = {'size': 50,
                            'color': 'black',
                            'multi': 'html'}

        node['label'] = f'<i>{node["id"]}</i>'

        #import pdb
        #pdb.set_trace()


    G_act = nx.from_pandas_edgelist(graph_df[graph_df.weight > 0], source='source', target='target', edge_attr='weight', create_using=nx.DiGraph())
    G_inact = nx.from_pandas_edgelist(graph_df[graph_df.weight < 0], source='source', target='target', edge_attr='weight', create_using=nx.DiGraph())


    return net, G_act, G_inact

def get_de_subset(timepoint, chamber):
    if timepoint == 'e85':
        if chamber == 'IFT':
            return ['Xist', 'Uty', 'Tsix', 'Tnnc1', 'Igfbp5', 'Pde4d', 'Hsp90aa1', 'Rps29', 'Rps28', 'Myocd', 'Hsp90ab1', 'Csrp3', 'Rpl39', 'Tmsb10', 'Rpl37', 'Rpl37a', 'Rps27', 'Eif2s3y', 'Mir99ahg', 'Rps2', 'Rpl41', 'Prkg1', 'Myl7', 'Ttn', 'Cacna2d2', 'Rps26', 'Tmsb4x', 'Ctsl', 'Podxl', 'Magi1', 'Rps3a1', 'Cacna1c', 'Myl3', 'Slc1a3', 'Nfia', 'Vcan', 'Rps21', 'Tenm3', 'Alcam', 'Slit3', 'Kdm5d', 'Rpl38', 'Nnat', 'Rpl29', 'Ddx3y', 'Fbxo32', 'Bsg', 'Myh6', 'Angpt1', 'Vim', 'Stmn1', '5430431A17Rik', 'Scd2', 'Gm10076', 'Rpl36', 'Epha7', 'Mgarp', 'Itpr1', 'Uba52', 'Malat1', 'Prickle1', 'Rrad', 'Ldb3', 'Armh4', 'Flrt3', 'Rpl35', 'Rpl36a', 'Tenm2', 'Wnt2', 'Nhs', 'mt-Atp6', 'Arid3b', 'Lbh', 'Ncl', 'Sparc', 'Hspd1', 'Gm47283', 'mt-Nd3', 'Meis2', 'Maged2', 'Rcsd1', 'Prtg', 'Cald1', 'Tenm4', 'Qk', 'Ddr1', 'Rpl34', 'Snrpf', 'Man1a', 'Setbp1', 'Bnip3', 'Tbx2', 'Nrp1', 'Cited2', 'Sox6', 'Crip2', 'Ppp1r9a', 'Myom1', 'Dchs2', 'Smad6', 'Kitl', 'Casz1', 'Mat2a', 'Snrpg', 'Tnni1', 'Ccbe1', 'Akt3', 'Ndn', 'Nap1l1', 'Ddx3x', 'Igf2', 'Atpif1', 'Zfpm2', 'Dach1', 'Fhod3', 'Igf1r', 'Syne2', 'Mybpc3', 'Hdac9', 'Ctsd', 'Rpa3', 'Aldoa', 'Cox7c', 'Col4a5', 'Vcam1', 'Atp5e', 'Mical2', 'Lclat1', 'Fkbp1a', 'Sec62', 'Nav2', 'Frmd4b', 'P3h2', 'Nebl', 'Fbn2', 'Jag1', 'Pdgfra', 'Reck', 'Celf2', 'Foxp1', 'Rtn4', 'Epb41l2', 'Pcdh7', 'Mfap2', 'Fbxl20', 'Cdk6', 'Crabp2', 'Tpm1', 'Atp8a1', 'Ppp1r14c', 'Npm1', 'Tomm7', 'Rps12', 'Vsnl1', 'Cox17', 'Flrt2', 'Robo2', 'Tusc3', 'Mab21l2', 'Mmp14', 'Edil3', 'Esrrg', 'Cdv3', 'Afdn', 'Kdm5b', 'Tmtc2', 'Rcn3', 'Kcnip1', 'Rps23', 'Ptn', 'Actc1', 'Rps20', 'mt-Nd4l', 'Gm49708', 'Gm26917', 'Rpl32', 'Ctnna3', 'Myh7', 'Col5a1', 'Hand2os1', 'Gata4', 'Ccng2', 'Eprs', '2610307P16Rik', 'Purb', 'Smyd1', 'Hsph1', 'Fbxl22', 'Nexn', 'Svil', 'Fau', 'Lrrfip1', 'Zbtb20', 'Atp5md', 'Pdlim5', 'Jarid2', 'S100a10', 'Mif', 'Eif3a', 'Map3k20']
        if chamber == 'V':
            return ['Xist', 'Tnnc1', 'Myl3', 'Actn2', 'Myom1', 'Myh7', 'Ttn', 'Myl7', 'Fbxl22', 'Csrp3', 'Mef2c', 'Hspb1', 'Robo2', 'Tenm3', 'Nrp1', 'Col2a1', 'Sox4', 'Asb2', 'Map3k20', 'Uty', 'Myl2', 'Cacna1c', 'Vim', 'Igfbp5', 'Myh6', 'Plxna4', 'Rrad', 'Ckb', 'Myo18b', 'Cald1', 'Pam', 'Mybpc3', 'Des', 'Cdkn1c', 'Basp1', 'Tpm1', 'Slc25a4', 'Kif26b', 'Scd2', 'Trdn', 'Hspb2', 'Gpc3', 'Smpx', 'Synpo2l', 'Zfpm2', 'Erc2', 'Fbxo32', 'Trim55', 'Marcks', 'Sh3bgr', 'Nnat', 'Tanc2', 'Xirp1', 'Ccbe1', '2610307P16Rik', 'Igf2', 'Wnt2', 'Rtn4', 'Armh4', 'Foxp1', 'Lbh', 'Igf2r', 'Myl4', 'Tmsb10', 'Fam49a', 'Pgam2', '5430431A17Rik', 'Auts2', 'Lama4', 'Tenm4', 'Thbs4', 'Magi1', 'Hspb7', 'Palld', 'Rps2', 'Esrrg', 'Tuba1a', 'Syne2', 'Fkbp1a', 'Cryab', 'Slc22a17', 'Eno3', 'Gata4', 'Chchd10', 'Crip2', 'Tceal9', 'Nfia', 'Ctsl', 'Ctnna3', 'Stmn1', 'Pcdh7', 'Pde4b', 'Cacna2d2', 'Tbx20', 'Bzw2', 'Mest', 'Flna', 'Tmem163', 'Angpt1', 'Tsix', 'Hcfc1r1', 'Hsp90b1', 'Ldb3', 'Kdm5b', 'Kitl', 'Obscn', 'Col18a1', 'Tnni1', 'Ak1', 'Ddr1', 'Hrc', 'Maged1', 'Lclat1', 'Man1a', 'Rps27', 'Smad6', 'Rpsa', 'Flnc', 'Malat1', 'Bbx', 'Actc1', 'Sorbs2', 'Csrp1', 'Hapln1', 'Pde4dip', 'Mmp14', 'Hs6st2', 'Rps29', 'Cacna1d', 'Rps28', 'Rpl39', 'Cobll1', 'Smyd1', 'Calr', 'Wwc2', 'Slc6a6', 'Dpysl3', 'Tmsb4x', 'Sptbn1', 'Abra', 'Mab21l2', 'Itga6', 'Atcayos', 'Frmd4b', 'Filip1', 'Tmem56', 'Sema3a', 'Gm49708', 'Mgarp', 'Ano4', 'Eif2s3y', 'Afdn', 'Mdk', 'Zbtb20', 'Rps26', 'Fabp3', 'Ephb1', 'Gulp1', 'Igf1r', 'Meis2', 'Tenm2', 'Resf1', 'Rpl41', 'Reep1', 'Lrrc10', 'Atp8a1', 'Edil3', 'Fsd2', 'Atpif1', 'Acta1', 'Tmod1', 'Nkx2-5', 'Slc1a3', 'Baz2b', 'Gm10076', 'Pik3r3', 'Ppp2r5a', 'Cemip2', 'Mical2', 'Smarcc1', 'Rpl7', 'Crybg3', 'Slc2a1', 'Kcnh7', 'Col5a1', 'Slc16a3', 'Myl9', 'Itpr1', 'Ndn', 'Rpl38', 'Mir99ahg', 'Bsg', 'Tcap', 'Fhl2', 'Ppp2r3a', 'Ddx3y', 'Ccng2', 'Pdgfra', 'Thsd4', 'Rpl13']
        
    else:
        if chamber == 'A':
            return ['Hba-x', 'Hbb-y', 'Hbb-bh1', 'mt-Nd1', 'Jarid2', 'mt-Co3', 'mt-Atp6', 'Xist', 'mt-Nd4l', 'mt-Nd3', 'mt-Nd2', 'Tnnc1', 'mt-Nd4', 'Hba-a1', 'mt-Cytb', 'mt-Co2', 'Myl3', 'Rpl29', 'Igf2r', 'AY036118', 'mt-Nd5', 'Myl7', 'Myl2', 'Grb10', 'Rpl27', 'Csrp3', 'Bsg', 'Tsix', 'Myom1', 'P4ha1', 'Gm42418', 'Atp5g1', 'Cycs', 'Ctnna3', 'Mef2c', 'Magi1', 'Fdps', 'Peg3', 'Sorbs2', 'Igfbp5', 'Ctsl', 'Cmss1', 'Slc2a1', 'Auts2', 'Hspb1', 'Myh7', 'Igf1r', 'Hk2', 'Ttn', 'Mybpc3', 'Mir99ahg', 'Slc25a5', 'Gm10076', 'Igf2', 'Slc2a3', 'Ppp1r14c', 'Hspa8', 'H2afz', 'Sfrp1', 'Foxp1', 'Plod2', 'Myh6', 'Actn2', 'Hs6st2', 'Celf2', 'Ccnd2', 'mt-Atp8', 'Uty', 'Fam162a', 'Sptbn1', 'Tnnt2', 'Ero1l', 'Dpp6', 'Cox8a', 'Crip2', 'Igfbp2', 'Gpc3', 'Atp5j', 'Rps12', 'Fras1', 'Pgk1', 'Vcan', 'Clcn3', '5430431A17Rik', 'Myo18b', 'Loxl2', 'Lpp', 'Rpl41', 'Nrp1', 'Colec12', 'Rplp1', 'Tmtc2', 'Wls', 'Actg1', 'Tnni1', 'Camk1d', 'Atp5g3', 'Cox6c', 'Rpl35', 'Cox5b', 'Smad6', 'Cald1', 'Actb', 'Tenm4', 'Pgam2', 'Cryab', 'Slc25a4', 'Tceal9', 'Csrp2', 'Cox5a', 'Nos1ap', 'Trim55', 'Plxna4', 'Fbxl22', 'Zfpm2', 'Nfia', 'Uqcrb', 'Atp5j2', 'Ldb3', 'Trdn', 'Gm10260', 'Cacna1c', 'Hsp90aa1', 'Atp5a1', 'Myl4', 'Myh10', 'Myl1', 'Tnni3', 'Rps23', 'Soat1', 'Ppia', 'Map3k20', 'Cyc1', 'Rps11', 'Chchd10', 'Bnip3', 'Hdac9', 'Filip1', 'mt-Co1', 'Zbtb20', 'Sh3bgr', 'Rrad', 'Nppa', 'Usp9x', 'Fkbp1a', 'Ppp1r9a', 'Asb2', 'Elovl6', 'Hspb2', 'Palld', 'Fam49a', 'Smyd1', 'Rps13', 'Hmgcs1', 'Esrrg', 'Edil3', 'Ndufa4', 'Uqcrh', 'Rps26', 'Idi1', 'Slc16a3', 'Tmsb10', 'Col2a1', 'Atp5b', 'Synpo2l', 'Tubb4b', 'Dhcr24', 'Obscn', 'Homer2', 'Flrt2', 'Bbx', 'Atp5o', 'Flnc', 'Prkca', 'Cpeb2', 'Clic4', 'Foxp2', 'Ndufa5', 'Malat1', 'Acaca', 'Peg10', 'Rplp2', 'Kdm5b', 'Smad7', 'Rps15a', 'Gm49708', 'Rbms1', 'P4ha2', 'Des', 'Hspb7', 'Gata4', 'Actc1', 'Xirp1', 'Rps27a', 'Rian', 'Atp5h', 'Epha7', 'Rpl17', 'Smpx', 'Akt3']
        if chamber == 'V':
            return ['Hba-x', 'Hbb-y', 'mt-Nd1', 'mt-Atp6', 'Myl3', 'Tnnc1', 'Hbb-bh1', 'mt-Nd3', 'mt-Nd2', 'mt-Nd4', 'mt-Co3', 'mt-Cytb', 'Jarid2', 'Trdn', 'Igf2r', 'mt-Co2', 'mt-Nd4l', 'Myh7', 'Magi1', 'Hspa8', 'Slc25a5', 'Myl7', 'Rpl29', 'Nfia', 'Slc8a1', 'Atp5g1', 'Nppa', 'Des', 'Actn2', 'Mef2c', 'Cycs', 'mt-Nd5', 'Smpx', 'Rbpms', 'Slc2a1', 'Ank3', 'Armh4', 'Myl2', 'Cryab', 'Cacna1d', 'Atp5b', 'Tenm3', 'Edil3', 'Lpp', 'Kcnh7', 'Hsp90aa1', 'Ccnd2', 'Hspd1', 'Msi2', 'Hba-a1', 'Casq1', 'Tbx20', 'Zbtb20', 'Tanc2', 'Atp5g3', 'Cox5b', 'mt-Atp8', 'Tpm1', 'Sptbn1', 'Atp5o', 'Rbm20', 'Cox8a', 'Cdkn1c', 'Esrrg', 'Robo2', 'Chchd10', 'Atp5h', 'Ndufa4', 'Atp5a1', 'Rps23', 'Mest', 'Ryr2', 'Rbfox2', 'Palld', 'Sh3bgr', 'Grb10', 'Syne2', 'Rpl27', 'Rabgap1l', 'Airn', 'AY036118', 'Xist', 'Actg1', 'Erbb4', 'Fbxo32', 'Cdh2', 'Wnt2', 'Marcks', 'Atp5j', 'Myom1', 'Foxp1', 'Afdn', 'Ctsl', 'Eno3', 'Colec12', 'Map3k20', 'Pdlim5', 'Prkg1', 'Malat1', 'Pard3', 'Ncl', 'Rplp1', 'Acta1', 'Gm49708', 'Arhgap6', 'Kdm5b', 'Sox4', 'Grip1', 'Lclat1', 'Zfpm2', 'Myl4', 'Celf2', 'Atp5j2', 'Igf1r', 'Obscn', 'Uqcrh', 'Rps8', 'Lama4', 'Col2a1', 'Pgam2', 'Cox6c', 'Gpc3', '2610307P16Rik', 'Cyc1', 'Hk2', 'Itga6', 'Ndufa5', 'Gm42418', 'Igfbp5', 'Slc25a4', 'Nrp1', 'Mical2', 'Cox5a', 'Rpl41', 'Slc6a6', 'Fkbp1a', 'Bnip3', 'Rps13', 'Qk', 'Hcfc1r1', 'Etfb', 'Hspb7', 'Ppia', 'Actb', 'Rere', 'Gata4', 'Tubb4b', 'Synpo2l', 'Mir99ahg', 'Cdh13', 'Uqcrb', 'Idh3a', 'Rps11', 'Rps15a', 'Pam', 'Hspb2', 'Fbxl22', 'Pfkp', 'Nasp', 'Rps12', 'Pde4b', 'Mitf', 'Rasgef1b', 'Myl1', 'Nos1ap', 'Rps27a', 'Cox6a2', 'Ctnna1', 'Rpl23', 'St6galnac3', 'mt-Co1', 'Gm10076', 'Tcap', 'Creb3l2', 'Pde3a', 'Tmtc2', 'Rps26', 'Mid2', 'Smad6', 'Hspb1', 'Tgfb2', 'Plod2', 'Rpl11', 'Rpl27a', 'Wnk1', 'Spon1', 'Sipa1l2', 'Rtn4', 'Pakap.1', 'Kdm7a', 'Mgarp', 'Meg3', 'Vim', 'Prkce', 'Cmss1', 'Id3', 'Rspo3', 'Wwc2', 'Akap9', 'Atp5d']

def main():
    wt_e85_ift = pd.read_csv('./data/tmp_e85-wt.csv')
    ko_e85_ift = pd.read_csv('./data/tmp_e85-ko.csv')
    
    wt_e85_v = pd.read_csv('./data/tmp_e85-wt-V.csv')
    ko_e85_v = pd.read_csv('./data/tmp_e85-ko-V.csv')
    
    wt_e9_a = pd.read_csv('./data/tmp_e9-wt-A.csv')
    ko_e9_a = pd.read_csv('./data/tmp_e9-ko-A.csv')
    
    wt_e9_v = pd.read_csv('./data/tmp_e9-wt-V.csv')
    ko_e9_v = pd.read_csv('./data/tmp_e9-ko-V.csv')
    
    # E85 IFT
    visualize_sub(wt_e85_ift, get_de_subset('e85', 'IFT'), 'wt_e85_ift.html')
    visualize_sub(ko_e85_ift, get_de_subset('e85', 'IFT'), 'ko_e85_ift.html')
    
    # E85 V
    visualize_sub(wt_e85_v, get_de_subset('e85', 'V'), 'wt_e85_v.html')
    visualize_sub(ko_e85_v, get_de_subset('e85', 'V'), 'ko_e85_v.html')
    
    # E9 A
    visualize_sub(wt_e9_a, get_de_subset('e9', 'A'), 'wt_e9_A.html')
    visualize_sub(ko_e9_a, get_de_subset('e9', 'A'), 'ko_e9_A.html')
    
    # E9 V
    visualize_sub(wt_e9_v, get_de_subset('e9', 'V'), 'wt_e9_V.html')
    visualize_sub(ko_e9_v, get_de_subset('e9', 'V'), 'ko_e9_V.html')
    


if __name__ == '__main__':
    main()
