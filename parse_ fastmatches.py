#to get the match results, in directory with fm result json files:
import os
import json
results=[i for i in  os.listdir() if i[-4:]=='json']
matches=[]
for f in results:
  df=open(f,'r')
  res=json.load(df)
  for e in res:#each one of these is a dictionary with scaff,motif,nmatches,matched pos
      x=str(e).split(',')
      nmatches=x[2].split(':')[1].strip()
      if int(nmatches)>0:
          matches.append(e)
          # for key in e.keys():
          #     print(str(key) + ': '+str(e[key]))
          # print('\n')
  df.close()

out_file = open("matches.json", "w")
json.dump(matches, out_file, indent = 6)
out_file.close()


'''
In [10]: len(matches)
Out[10]: 14412

you know what, this along with the degeneracy of matched positions leads me to
think maybe some of my 2500 motifs are the same
EDIT FILTER MOTIF SCRIPT TO ENSURE NONREDUNDANCY

for now gonna write these to a file to easily load up later
ofile=open('allmatches.txt','w')
ofile.write(str(matches))
ofile.close()
'''


df=open("matches.json",'r')
matches=json.load(df)
df.close()

#okay make a new directory that will be used for standard matching
os.mkdir('rosettamatch')
os.chdir('rosettamatch')
os.mkdir('cluster_out')
#lets generate motif csts and put them in this directory
#along with position files
#for which i am only using fast matched positions here
#append cst and pos file names to matches directories
paramspath='/wynton/home/kortemme/cgalvin/newtargets/nps/Inputs/Rosetta_Inputs/nps.params' #####
shellfile_suffix='.sh'
#
count=1
indices=[i for i in range(0,len(matches),100)]
indices.append(len(matches))
for i,e in enumerate(indices[:-1]):
    sfn='motif2cst_'+str(count)+shellfile_suffix
    sf=open(sfn,'w')
    sf.write('#!/bin/bash\n')
    sf.write('source ~/anaconda3/etc/profile.d/conda.sh\n')
    sf.write('conda activate pyr37\n')
    for d in matches[e:indices[i+1]]:
        motif_path=d['motif']
        scaffold_name=d['scaffold'].split('/')[-1].strip('.pdb.gz')
        motif_name=motif_path.split('/')[-1].strip('.pdb')
        cst_name=scaffold_name+'_match_'+motif_name+'.cst'
        positions=d['fast_matched_positions']
        pos_file_name=scaffold_name+'_match_'+motif_name+'.pos'
        of=open(pos_file_name,'w')
        pos_string=''
        for l in positions:
            for i in l:
                pos_string+=str(i)+'  '
        of.write(pos_string)
        of.close()
        s='time python ~/tools/motif_to_cst.py '+motif_path+' '+paramspath+' '+cst_name
        d['cst']=cst_name
        d['pos']=pos_file_name
        sf.write(s+'\n')
    count+=1
    sf.close()
    os.system('chmod ugo+x '+sfn)
    os.system('qsub -cwd -l mem_free=2G -e cluster_out -o cluster_out '+sfn)

out_file = open("matches.json", "w")
json.dump(matches, out_file, indent = 6)
out_file.close()
'''
count=1
indices=[i for i in range(0,len(matches),100)]
indices.append(len(matches))
shellfile_suffix='.sh'
for i,e in enumerate(indices[:-1]):
    sfn='motif2cst_'+str(count)+shellfile_suffix
    os.system('chmod ugo+x '+sfn)
    os.system('qsub -cwd -l mem_free=2G -e cluster_out -o cluster_out '+sfn)
    count+=1

#to just do from terminal and not submit motif2cst jobs
for d in matches:
    motif_path=d['motif']
    scaffold_name=d['scaffold'].split('/')[-1].strip('.pdb.gz')
    motif_name=motif_path.split('/')[-1].strip('.pdb')
    cst_name=scaffold_name+'_match_'+motif_name+'.cst'
    positions=d['fast_matched_positions']
    pos_file_name=scaffold_name+'_match_'+motif_name+'.pos'
    of=open(pos_file_name,'w')
    pos_string=''
    for l in positions:
        for i in l:
            pos_string+=str(i)+'  '
    of.write(pos_string)
    of.close()
    s='time python ~/tools/motif_to_cst.py '+motif_path+' '+paramspath+' '+cst_name
    os.system(s)
    d['cst']=cst_name
    d['pos']=pos_file_name
'''
import os
import json
df=open("matches.json",'r')
matches=json.load(df)
df.close()
os.chdir('rosettamatch')
##################################
#okay now I can set up matcher job
#now to run the standard matcher
##################################
#parameters to specify
##################################
paramspath='/wynton/home/kortemme/cgalvin/newtargets/nps/Inputs/Rosetta_Inputs/nps.params' #####
scaffold_path_prefix='/wynton/home/kortemme/cgalvin/lucs_scaffolds/'
# motif_path_prefix='/wynton/home/kortemme/cgalvin/bsff_motifs/nps_2_motif_pdbs/clean_motifs/nps_filtered_motifs/'
cst_pos_path_prefix='/wynton/home/kortemme/cgalvin/lucs_scaffolds/nps_2lv8_fm_output/rosettamatch/'
lig_name='nps'
rmatch_output_dir='nps_standard_matches'
os.mkdir(rmatch_output_dir)
shellfile_suffix='_nps_standard_match.sh'
##################################
#submit
##################################
count=1
indices=[i for i in range(0,len(matches),100)]
indices.append(len(matches))
for i,e in enumerate(indices[:-1]):
    sfn='_'+str(count)+shellfile_suffix
    of=open(sfn,'w')
    of.write('#!/bin/bash\n')
    for d in matches[e:indices[i+1]]:
        scaffold_path=scaffold_path_prefix+d['scaffold']
        motif_path=d['motif']
        pos_path=cst_pos_path_prefix+d['pos']
        cst_path=cst_pos_path_prefix+d['cst']
        matcher_cmd = ['~/main/source/bin/match.linuxgccrelease',
                       '-s', scaffold_path,
                       '-extra_res_fa', paramspath,
                       '-match:geometric_constraint_file', cst_path,
                       '-match::scaffold_active_site_residues', pos_path,
                       '-match:output_format', 'PDB',
                       '-match:match_grouper', 'SameSequenceGrouper',
                       '-match:consolidate_matches',
                       '-match:output_matches_per_group', '1',
                       '-use_input_sc',
                       # '-in:ignore_unrecognized_res',
                       '-ex1', '-ex2','-extrachi_cutoff 0',
                       '-enumerate_ligand_rotamers', 'false',
                       '-match::lig_name', lig_name,
                       '-out:path:all',rmatch_output_dir]
        cmd=' '.join(matcher_cmd)
        of.write(cmd+'\n')
    count+=1
    of.close()
    os.system('chmod ugo+x '+sfn)
    os.system('qsub -cwd -l mem_free=2G '+sfn)




# some jobs crash, how many?

f=open('ROSETTA_CRASH.log','r' )
lines=[line for line in f.readlines()]
f.close()
count=0
for i in lines:
    if i=='[START_CRASH_REPORT]\n':
        count+=1

'''
In [2]: count
Out[2]: 4015

okay 80% success isnt bad, ill proceed to submit all 14k

'[START_MESSAGE]\n',
 'Cannot create normalized xyzVector from vector of length() zero.\n',
'''

standardmatches=[i for i in os.listdir() if i[:2]=='UM']
'''
In [6]: len(standardmatches)
Out[6]: 16
'''

scaffold_path_prefix='/wynton/home/kortemme/cgalvin/lucs_scaffolds/'
motif_path_prefix='/wynton/home/kortemme/cgalvin/bsff_motifs/nps_2_motif_pdbs/clean_motifs/nps_filtered_motifs/'

os.mkdir('nps_matches')
for match in standardmatches:
    motif=motif_path_prefix+match.split('_match_')[1][:-6]+'.pdb'
    # scaffold=scaffold_path_prefix+match.split('_match_')[0].split('_')[4]+'_'+(match.split('_match_')[0].split('_')[5])+'.pdb.gz'
    os.system('cp '+motif+' nps_matches/'+motif.split('/')[-1])
    # os.system('cp '+scaffold+' nps_matches/'+scaffold.split('/')[-1])
    os.system('cp '+match+' nps_matches/'+match)
#scaffold needs number of subdir in 2lv8 dir

'''
Executive: RMS =    0.769 (64 to 64 atoms)
Executive: RMS =    0.018 (16 to 16 atoms)
'''











################################################
################################################
################################################
################################################
#SUBMITTING STANDARD MATCHER JOB WITHOUT FASTMATCH
################################################
################################################
################################################
################################################

#in filtered motif directory:
import os
motif_path='/wynton/home/kortemme/cgalvin/bsff_motifs/3ng_2_motif_pdbs/clean_motifs/3ng_filtered_motifs'
motifs=[i for i in os.listdir(motif_path) if i[-3:]=='pdb']

os.mkdir('rosettamatch')
os.chdir('rosettamatch')
os.mkdir('cluster_out')
#lets generate motif csts and put them in this directory
paramspath='/wynton/home/kortemme/cgalvin/newtargets/3ng/Inputs/Rosetta_Inputs/3ng.params' #####
shellfile_suffix='.sh'
count=1
for motif in motifs:
    sfn='motif2cst_'+str(count)+shellfile_suffix
    sf=open(sfn,'w')
    sf.write('#!/bin/bash\n')
    sf.write('source ~/anaconda3/etc/profile.d/conda.sh\n')
    sf.write('conda activate pyr37\n')
    cst_name=motif.split('.')[0]+'.cst'
    s='time python ~/tools/motif_to_cst.py '+os.path.join(motif_path,motif)+' '+paramspath+' '+cst_name
    sf.write(s+'\n')
    count+=1
    sf.close()
    os.system('chmod ugo+x '+sfn)
    os.system('qsub -cwd -e cluster_out -o cluster_out '+sfn)

'''
#now lets generate pos files and also put them in this directory
#get scaffold pos files in matches dir
import os
scaffold_path='/wynton/home/kortemme/cgalvin/lucs_scaffolds/filtered'
scaffolds=[i for i in os.listdir(scaffold_path) if i[-3:]=='.gz' or i[-3:]=='pdb']

In [7]: len(scaffolds)
Out[7]: 1506

#putting pos in filtered scaffolds dir so i should only have to do it this
#one time now ......

######REDO IF MORE FILTERED SCAFFOLDS ADDED!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from pyrosetta import *
init('-ignore_unrecognized_res')
for scaffold in scaffolds:
    posname=scaffold.split('.')[0]+'.pos'
    p=pose_from_pdb(os.path.join(scaffold_path,scaffold))
    nres=p.total_residue()
    of=open(os.path.join(scaffold_path,posname),'w')
    for i in range(nres):
        of.write(str(i+1)+'  ')
    of.close()
'''

#now to actually submit matcher jobs
#each job is one motif with all scaffolds
#in rosettamatch directory:
import os
paramspath='/wynton/home/kortemme/cgalvin/newtargets/3ng/Inputs/Rosetta_Inputs/3ng.params' #####
scaffold_path='/wynton/home/kortemme/cgalvin/lucs_scaffolds/filtered'
scaffolds=[i for i in os.listdir(scaffold_path) if i[-3:]=='.gz' or i[-3:]=='pdb']
motif_path='/wynton/home/kortemme/cgalvin/bsff_motifs/3ng_2_motif_pdbs/clean_motifs/3ng_filtered_motifs'
motifs=[i for i in os.listdir(motif_path) if i[-3:]=='pdb']
lig_name='3ng'
shellfile_suffix='_3ng_match.sh'
##################################
#submit
##################################
count=1
for motif in motifs:
    sfn='_'+str(count)+shellfile_suffix
    of=open(sfn,'w')
    of.write('#!/bin/bash\n')
    cst_path=motif.split('.')[0]+'.cst'
    for scaffold in scaffolds:
        scaffold_input=os.path.join(scaffold_path,scaffold)
        pos_name=scaffold.split('.')[0]+'.pos'
        pos_path=os.path.join(scaffold_path,pos_name)
        matcher_cmd = ['~/main/source/bin/match.linuxgccrelease',
                       '-s', scaffold_input,
                       '-extra_res_fa', paramspath,
                       '-match:geometric_constraint_file', cst_path,
                       '-match::scaffold_active_site_residues', pos_path,
                       '-match:output_format', 'PDB',
                       '-match:match_grouper', 'SameSequenceGrouper',
                       '-match:consolidate_matches',
                       '-match:output_matches_per_group', '1',
                       '-use_input_sc',
                       '-ex1', '-ex2','-extrachi_cutoff 0',
                       '-enumerate_ligand_rotamers', 'false',
                       '-match::lig_name', lig_name]
        cmd=' '.join(matcher_cmd)
        of.write(cmd+'\n')
    of.write('\nqstat -j "$JOB_ID"')
    count+=1
    of.close()
    os.system('chmod ugo+x '+sfn)
    os.system('qsub -cwd -l mem_free=2G -o cluster_out -e cluster_out '+sfn)



#########
#get matches and consolidate to one directory
import os
standardmatches=[i for i in os.listdir() if i[:2]=='UM']
'''
In [3]: len(standardmatches)
Out[3]: 26
'''
resultdir='nps_matches2'
os.mkdir(resultdir)
motif_path='/wynton/home/kortemme/cgalvin/bsff_motifs/nps_2_motif_pdbs/clean_motifs/nps_filtered_motifs/'
for match in standardmatches:
    motifres=match.split('_')[6:9]
    motif=motif_path+'_'.join(motifres)+'.pdb'
    # scaffold=scaffold_path_prefix+match.split('_match_')[0].split('_')[4]+'_'+(match.split('_match_')[0].split('_')[5])+'.pdb.gz'
    os.system('cp '+motif+' '+resultdir+'/'+motif.split('/')[-1])
    # os.system('cp '+scaffold+' nps_matches/'+scaffold.split('/')[-1])
    os.system('cp '+match+' '+resultdir+'/'+match)
#scaffold needs number of subdir in 2lv8 dir

'''
nps WOW SOMETHING FUNKY HAPPENED WITH THE CLUSTER AND IT CORRUPTED MY FILES ONLY 13 NOW
PERHAPS I CAN RECOVER THEM BY SIMPLY RESUBMITTING THE JOBS USING THE TITLE OF THE FILES
TO ID APPROPRIATE SCAFFOLD AND MOTIF INOUT
RESUBMITTING THESE CORRUPTED JOBS:

DEX, A8S, 3NG ALL HAVE WEIRD CORRUPT MATCH PROBLEM, RESUBMITTING
            3ng doesnt even have match files though....

import os
standardmatches=[i for i in os.listdir() if i[:2]=='UM']
corrupt=[]
for i in standardmatches:
    f=open(i,'r')
    lines=[line for line in f.readlines()]
    f.close()
    if len(lines)==0:
        corrupt.append(i)
os.mkdir('matches')
for i in standardmatches:
    f=open(i,'r')
    lines=[line for line in f.readlines()]
    f.close()
    if len(lines)==0:
        corrupt.append(i)
    else:
        if i not in os.listdir('matches'):
            os.system('cp '+i+' matches/'+i)

for i in corrupt:
    os.system('mv '+i+' corrupt_'+i)

for i in os.listdir():
    if i[:7]=='corrupt':
        corrupt.append(i)

scaffold_path_prefix='/wynton/home/kortemme/cgalvin/lucs_scaffolds/filtered/'
motif_path_prefix='/wynton/home/kortemme/cgalvin/bsff_motifs/a8s_2_motif_pdbs/clean_motifs/a8s_filtered_motifs/'
paramspath='/wynton/home/kortemme/cgalvin/newtargets/a8s/Inputs/Rosetta_Inputs/a8s.params' #####
lig_name='a8s'
cst_prefix='/wynton/home/kortemme/cgalvin/bsff_motifs/a8s_2_motif_pdbs/clean_motifs/a8s_filtered_motifs/rosettamatch/'
shellfile_suffix='_a8s_match.sh'
count=1
for i in corrupt:
    if i[:7]=='corrupt':
        scaffold=i.split('_')[5:7]
    else:
        scaffold=i.split('_')[4:6]
    # print(scaffold)
    if '_'.join(scaffold)[:5]=='clean':
        scaffold_path=os.path.join(scaffold_path_prefix,'_'.join(scaffold)+'.pdb')
    else:
        scaffold_path=os.path.join(scaffold_path_prefix,'_'.join(scaffold)+'.pdb.gz')
    # print(scaffold_path)
    if i[:7]=='corrupt':
        motif=i.split('_')[7:10]
    else:
        motif=i.split('_')[6:9]
    motif_path=os.path.join(motif_path_prefix,'_'.join(motif)+'.pdb')
    # print(motif_path)
    cst_path=os.path.join(cst_prefix,'_'.join(motif)+'.cst')
    pos_path=os.path.join(scaffold_path_prefix,'_'.join(scaffold)+'.pos')
    # print(cst_path)
    # print(pos_path)
    if os.path.exists(scaffold_path):
        # print('scaff')
        if os.path.exists(paramspath):
            # print('params')
            if os.path.exists(cst_path):
                # print('cst')
                if os.path.exists(pos_path):
                    # print('pos')
                    sfn='_redo'+str(count)+shellfile_suffix
                    of=open(sfn,'w')
                    of.write('#!/bin/bash\n')
                    matcher_cmd = ['~/main/source/bin/match.linuxgccrelease',
                                   '-s', scaffold_path,
                                   '-extra_res_fa', paramspath,
                                   '-match:geometric_constraint_file', cst_path,
                                   '-match::scaffold_active_site_residues', pos_path,
                                   '-match:output_format', 'PDB',
                                   '-match:match_grouper', 'SameSequenceGrouper',
                                   '-match:consolidate_matches',
                                   '-match:output_matches_per_group', '1',
                                   '-use_input_sc',
                                   '-ex1', '-ex2','-extrachi_cutoff 0',
                                   '-enumerate_ligand_rotamers', 'false',
                                   '-match::lig_name', lig_name]
                    cmd=' '.join(matcher_cmd)
                    of.write(cmd+'\n')
                    of.write('\nqstat -j "$JOB_ID"')
                    count+=1
                    of.close()
                    os.system('chmod ugo+x '+sfn)
                    os.system('qsub -cwd -l mem_free=2G -o cluster_out -e cluster_out '+sfn)


NOW MOVE THE REACQUIRED MATCHES TO MATCHES DIR


NEXT THING HERE WILL BE TO MAYBE TEST RELAXED CONSTRAINTS OR POTENTIALLY MORE
SCAFF/MOTIF DEPENDING ON HOW LONG THESE JOBS TAKE TO COMPLETE
'''


'''
MATCHING ERROR I KEEP SEEING
[START_MESSAGE]
Cannot create normalized xyzVector from vector of length() zero.

[END_MESSAGE]
[END_CRASH_REPORT]

FROM ROSETTACOMMONS:

Usually this means you have three colinear atoms for which the system is trying to compute a dihedral.

Look at your inputs and parameters to see if you've got three bonded colinear atoms, or bonded atoms that are defined with the same input coordinates the PDB file.  You can fuzz the coordinates a bit so they aren't perfectly colinear, and it should start working.

(Even if you are supposed to have a bunch of perfectly colinear atoms - fuzz them anyway, in the PDB and params if necessary.  A 1% error in colinearity is enough to make dihedral math work and below the resolution of the PDB file anyway).
'''



################################################
################################################
################################################
######################DESIGN####################
################################################
################################################


#submitting matches to design
import os
standardmatches=[i for i in os.listdir() if i[:2]=='UM']
#consolidate to a new folder 'matches'
# os.mkdir('matches')
corrupt=[]
count=0
for i in standardmatches:
    f=open(i,'r')
    lines=[line for line in f.readlines()]
    f.close()
    if len(lines)==0:
        corrupt.append(i)
    else:
        if i not in os.listdir('matches'):
            count+=1
            os.system('cp '+i+' matches/'+i)
print(str(count)+' proper new matches')
print(str(len(corrupt))+' new corrupt')

#enter this matches directory
#now get params for each and put in new directory with corrected pdb
import os
standardmatches=[i for i in os.listdir() if i[:2]=='UM']
#first I need to get gp params files and replace lig coords in pdb with mol2params output pdb
ligand_charge='-1'
ligname='nps'
cleandirname='design'
os.mkdir(cleandirname)
for initial_match in standardmatches:
    fname=initial_match.split('.')[0]
    f=open(initial_match,'r')
    lines=[]
    alines=[]
    for line in f.readlines():
        if line[0:6]=='HETATM':
            lines.append(line)
        elif line[0:4]=='ATOM':
            alines.append(line)
    f.close()
    newligpdb=open(fname+'_lig.pdb','w')
    for line in lines:
        newligpdb.write(line)
    newligpdb.close()
    os.system('obabel -i pdb '+fname+'_lig.pdb'+' -o mol2 -O '+fname+'_lig.mol2')
    os.system('obabel -i mol2 '+fname+'_lig.mol2'+' -o mol2 -O '+fname+'_ligH.mol2 -p 7.4')
    os.system('obabel -i mol2 '+fname+'_ligH.mol2'+' -o pdb -O '+fname+'_ligH.pdb')
    os.system('time antechamber -i '+fname+'_ligH.pdb'+' -fi pdb -o '+fname+'_ligH_bcc.mol2'+' -fo mol2 -c bcc -nc '+ligand_charge)
    os.system('~/desktop/rosetta/main/source/scripts/python/public/generic_potential/mol2genparams.py --nm '+ligname+' -s '+fname+'_ligH_bcc.mol2')
    newlig=open(fname+'_ligH_bcc_0001.pdb','r')
    newliglines=[line for line in newlig.readlines() if line[:6]=='HETATM']
    newpdb=open(os.path.join(cleandirname,fname+'_clean.pdb'),'w')
    for line in alines:
        newpdb.write(line)
    newpdb.write('TER\n')
    for line in newliglines:
        newpdb.write(line)
    newpdb.write('TER\n')
    newpdb.close()
    os.system('mv '+fname+'_ligH_bcc.params '+os.path.join(cleandirname,fname+'_ligH_bcc.params'))
    os.system('mv '+fname+'_ligH_bcc_0001.pdb '+os.path.join(cleandirname,fname+'_ligH_bcc_0001.pdb'))



import os
#have to fix the params files to make sure no redundant bonds defined :(
#again not entirely sure whats going on here or if i should change the 1 to 2s on these redundant lines
#as perhaps that defines as double bond? just leaving as 1 for now since my
#previous test of whether it made a difference showed that it didnt
#(ie got same exact score using different versions of params)
prms=[i for i in os.listdir() if i[-6:]=='params']
for prm in prms:
    f=open(prm,'r')
    lines=[line for line in f.readlines()]
    f.close()
    batoms=[]
    for ind,line in enumerate(lines):
        if line[:9]=='BOND_TYPE':
            a1=line[10:13].strip(' ')
            a2=line[15:18].strip(' ')
            batoms.append(((a1,a2),ind))
    toremove=[]
    for ix,e1 in enumerate(batoms):
        for e2 in batoms[ix+1:]:
            if e1!=e2:
                if e1[0]==e2[0]:
                    toremove.append(e2[1])
    toremove=sorted(toremove,reverse=True)
    for idx in toremove:
        lines.pop(idx)
    of=open(prm,'w')
    for line in lines:
        of.write(line)
    of.close()

#now get resfiles of new design inputs
#if im allowing lig to move, as cm does, then it doesnt make sense to
#constrain motif residues to native rotamer, so i am allowing them to repack
#in designs folder:
matches=[i for i in os.listdir() if i[:2]=='UM' and i[-9:]=='clean.pdb']
for match in matches:
    params = match.split('.')[0][:-6]+'_ligH_bcc.params'
    cmd='time python ~/desktop/match_resfile.py '+match+' '+params+' 5.0 9.0 genpot'
    print(cmd)
    os.system(cmd)


#now submit the cm jobs
import os
matches=[i for i in os.listdir() if i[:2]=='UM' and i[-9:]=='clean.pdb']
count=3
sfsuffix='nps_cm.sh'
for match in matches[2:]:
    sfn='_'+str(count)+sfsuffix
    sf=open(sfn,'w')
    params = match.split('.')[0][:-6]+'_ligH_bcc.params'
    resfile=match.split('.')[0]+'.resfile'
    outputdir='match_'+str(count)+'_designs'
    os.mkdir(outputdir)
    sf.write('#!/bin/bash')
    sf.write('\n')
    sf.write('tasks=(0\n')
    for x in range(1,100):
        sf.write('       '+str(x)+'\n')
    sf.write('       '+'100)')
    sf.write('\n~/main/source/bin/coupled_moves.default.linuxgccrelease -s '+match+' -mute protocols.backrub.BackrubMover -ex1 -ex2 -extrachi_cutoff 0 -nstruct 1 -score::weights beta_genpot -corrections::gen_potential -extra_res_fa '+params+' -resfile '+resfile+' -coupled_moves::mc_kt 2.4 -coupled_moves::boltzmann_kt 2.4 -coupled_moves::ntrials 1000 -coupled_moves::initial_repack false -coupled_moves::ligand_mode true -coupled_moves::ligand_weight 2 -coupled_moves::fix_backbone false -coupled_moves::bias_sampling true -coupled_moves::bump_check true -coupled_moves::exclude_nonclashing_positions true -coupled_moves::backbone_mover kic -coupled_moves::kic_perturber walking -nstruct 1 -gen_potential -score:set_weights hbond_bb_sc 2.0 -score:set_weights hbond_sc 2.0 -load_PDB_components False -out:path:all '+outputdir+' -out:prefix ${tasks[$SGE_TASK_ID]}_')
    sf.write('\nqstat -j "$JOB_ID"')
    sf.close()
    os.system('chmod ugo+x '+sfn)
    os.system('qsub -cwd -t 1-100 -l mem_free=2G -o cluster_out -e cluster_out '+sfn)
    count+=1

################################################
################################################
################################################
######################DESIGN###################
####################with constraints############
################################################
################################################
'''
get rmsd motifs and matches using
motif_to_match_rmsd.py

made xml for match analysis

import os
matches=[i for i in os.listdir() if i[:2]=='UM' and i[-9:]=='clean.pdb']

for match in matches:
    ms=match.split('_')
    mres=ms[:10]
    params='_'.join(mres)+'_ligH_bcc.params'
    if os.path.exists(params):
        s='~/desktop/rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease -s '+match+' -parser:protocol ~/desktop/analysis_scripts/match_analysis_genpot.xml -parser:view -in:file::extra_res_fa '+params+' -score::weights beta_genpot -corrections::gen_potential -load_PDB_components False -out:file:score_only matchscores.sc'
        os.system(s)
'''

HEADER                                            03-MAR-20   XXXX
EXPDTA    THEORETICAL MODEL
REMARK

        NEED TO REDO SAME BUT SIMPLY KEEP REMARK LINES IN PDB
        AND USE NEW SCRIPT TO GENERATE GP CST
        AND THEN SUBMIT WITH NEW XML AND COMMAND
        i also think i can expand design residues in resfile
        and test without exclude nonclashing too


#enter directory with matches that i wanna design
#now get params for each and put in new directory with corrected pdb
import os
standardmatches=[i for i in os.listdir() if i[:2]=='UM']
#first I need to get gp params files and replace lig coords in pdb with mol2params output pdb
ligand_charge='-1'
ligname='nps'
cleandirname='design'
os.mkdir(cleandirname)
for initial_match in standardmatches:
    fname=initial_match.split('.')[0]
    f=open(initial_match,'r')
    lines=[]
    alines=[]
    for line in f.readlines():
        if line[0:6]=='HETATM':
            lines.append(line)
        elif line[0:6]=='HEADER':
            alines.append(line)
        elif line[0:6]=='REMARK':
            alines.append(line)
        elif line[0:4]=='ATOM':
            alines.append(line)
    f.close()
    newligpdb=open(fname+'_lig.pdb','w')
    for line in lines:
        newligpdb.write(line)
    newligpdb.close()
    os.system('obabel -i pdb '+fname+'_lig.pdb'+' -o mol2 -O '+fname+'_lig.mol2')
    os.system('obabel -i mol2 '+fname+'_lig.mol2'+' -o mol2 -O '+fname+'_ligH.mol2 -p 7.4')
    os.system('obabel -i mol2 '+fname+'_ligH.mol2'+' -o pdb -O '+fname+'_ligH.pdb')
    os.system('time antechamber -i '+fname+'_ligH.pdb'+' -fi pdb -o '+fname+'_ligH_bcc.mol2'+' -fo mol2 -c bcc -nc '+ligand_charge)
    os.system('~/desktop/rosetta/main/source/scripts/python/public/generic_potential/mol2genparams.py --nm '+ligname+' -s '+fname+'_ligH_bcc.mol2')
    newlig=open(fname+'_ligH_bcc_0001.pdb','r')
    newliglines=[line for line in newlig.readlines() if line[:6]=='HETATM']
    newpdb=open(os.path.join(cleandirname,fname+'_clean.pdb'),'w')
    for line in alines:
        newpdb.write(line)
    newpdb.write('TER\n')
    for line in newliglines:
        newpdb.write(line)
    newpdb.write('TER\n')
    newpdb.close()
    os.system('mv '+fname+'_ligH_bcc.params '+os.path.join(cleandirname,fname+'_ligH_bcc.params'))
    os.system('mv '+fname+'_ligH_bcc_0001.pdb '+os.path.join(cleandirname,fname+'_ligH_bcc_0001.pdb'))



#enter design directory
import os
#have to fix the params files to make sure no redundant bonds defined :(
#again not entirely sure whats going on here or if i should change the 1 to 2s on these redundant lines
#as perhaps that defines as double bond? just leaving as 1 for now since my
#previous test of whether it made a difference showed that it didnt
#(ie got same exact score using different versions of params)
prms=[i for i in os.listdir() if i[-6:]=='params']
for prm in prms:
    f=open(prm,'r')
    lines=[line for line in f.readlines()]
    f.close()
    batoms=[]
    for ind,line in enumerate(lines):
        if line[:9]=='BOND_TYPE':
            a1=line[10:13].strip(' ')
            a2=line[15:18].strip(' ')
            batoms.append(((a1,a2),ind))
    toremove=[]
    for ix,e1 in enumerate(batoms):
        for e2 in batoms[ix+1:]:
            if e1!=e2:
                if e1[0]==e2[0]:
                    toremove.append(e2[1])
    toremove=sorted(toremove,reverse=True)
    for idx in toremove:
        lines.pop(idx)
    of=open(prm,'w')
    for line in lines:
        of.write(line)
    of.close()

#now get resfiles of new design inputs
#if im allowing lig to move, as cm does, then it doesnt make sense to
#constrain motif residues to native rotamer, so i am allowing them to repack
#in designs folder:
matches=[i for i in os.listdir() if i[:2]=='UM' and i[-9:]=='clean.pdb']
for match in matches:
    params = match.split('.')[0][:-6]+'_ligH_bcc.params'
    cmd='time python ~/desktop/analysis_scripts/match_resfile.py '+match+' '+params+' 5.0 10.0 genpot'
    print(cmd)
    os.system(cmd)

#now get new cst files that are compatible with genpot ligand
import os
matches=[i for i in os.listdir() if i[:2]=='UM' and i[-9:]=='clean.pdb']
for match in matches:
    ms=match.split('_')
    mres=ms[:10]
    params='_'.join(mres)+'_ligH_bcc.params'
    csss=match.split('.')[0]
    cstname=csss+'.cst'
    cmd='time python ~/desktop/motif_to_cst_genpot.py '+' '.join([match,params,cstname])
    print(cmd)
    os.system(cmd)


#now submit the cm jobs
import os
matches=[i for i in os.listdir() if i[:2]=='UM' and i[-9:]=='clean.pdb']
count=2
sfsuffix='nps_cm_cst.sh'
os.makedirs('cluster_out',exist_ok=True)
for match in matches[1:]:
    sfn='_'+str(count)+sfsuffix
    sf=open(sfn,'w')
    params = match.split('.')[0][:-6]+'_ligH_bcc.params'
    resfile=match.split('.')[0]+'.resfile'
    csss=match.split('.')[0]
    cstname=csss+'.cst'
    outputdir='match_'+str(count)+'_designs'
    os.mkdir(outputdir)
    sf.write('#!/bin/bash')
    sf.write('\n')
    sf.write('tasks=(0\n')
    for x in range(1,100):
        sf.write('       '+str(x)+'\n')
    sf.write('       '+'100)')
    sf.write('\n~/main/source/bin/rosetta_scripts.default.linuxgccrelease -s '+match+' -parser:protocol ~/coupled_moves.xml -parser:view -run:preserve_header -in:file::extra_res_fa '+params+' -resfile '+resfile+' -mute protocols.backrub.BackrubMover -ex1 -ex2 -extrachi_cutoff 0 -nstruct 1 -score::weights beta_genpot -corrections::gen_potential -coupled_moves::mc_kt 2.4 -coupled_moves::boltzmann_kt 2.4 -coupled_moves::ntrials 1000 -coupled_moves::initial_repack false -coupled_moves::ligand_mode true -coupled_moves::ligand_weight 2 -coupled_moves::fix_backbone false -coupled_moves::bias_sampling true -coupled_moves::bump_check true -coupled_moves::backbone_mover kic -coupled_moves::kic_perturber walking -coupled_moves::exclude_nonclashing_positions true -nstruct 1 -gen_potential -score:set_weights hbond_bb_sc 2.0 -score:set_weights hbond_sc 2.0 -load_PDB_components False -enzdes:cstfile '+cstname+' -out:path:all '+outputdir+' -out:suffix -out:prefix ${tasks[$SGE_TASK_ID]}')
    sf.write('\nqstat -j "$JOB_ID"')
    sf.close()
    os.system('chmod ugo+x '+sfn)
    os.system('qsub -cwd -t 1-100 -l mem_free=2G -o cluster_out -e cluster_out '+sfn)
    count+=1



################################################
################################################
################################################
######################DESIGN ANALYSIS###########
######################DESIGN ANALYSIS###########
######################DESIGN ANALYSIS###########
################################################
################################################
import os
design_dirs=[]
for i in os.listdir():
    if os.path.isdir(i)==True:
        if i[:5]=='match':
            design_dirs.append(i)

for designs_dir in design_dirs:
    os.chdir(designs_dir)
    os.system('time python ~/tools/pdblist.py pdblist.txt')
    os.chdir('..')
count=1
sfsuffix='_nps_cm_cst_analysis.sh'
for designs_dir in design_dirs:
    print(str(count))
    os.chdir(designs_dir)
    try:
        os.system('rm *.sh')
        pdbs=[i for i in os.listdir() if i[-3:]=='pdb']
        pdbex=pdbs[0]
    except:
        print(designs_dir)
        os.chdir('..')
        continue
    pdbexprcs=pdbex.split('_')[1:10]
    paramsn='UM_'+'_'.join(pdbexprcs)+'_ligH_bcc.params'
    print(paramsn)
    params='/wynton/home/kortemme/cgalvin/nps_cm_cst/'+paramsn
    if os.path.exists(params):
        print('good')
    else:
        print(designs_dir)
    outputdir='anl'
    os.makedirs(outputdir,exist_ok=True)
    scorefilename='nps_cm_cst_analysis_'+str(count)+'.sc'
    sfn='_'+str(count)+sfsuffix
    sf=open(sfn,'w')
    sf.write('#!/bin/bash')
    sf.write('\n')
    sf.write('tasks=(0\n')
    f=open('pdblist.txt','r')
    pdbs=[i.strip('\n') for i in f.readlines()]
    f.close()
    for x in pdbs[:-1]:
        sf.write('       '+str(x)+'\n')
    sf.write('       '+str(pdbs[-1])+')\n')
    sf.write('~/main/source/bin/rosetta_scripts.default.linuxgccrelease -nstruct 1 -jd2:ntrials 1 -s ${tasks[$SGE_TASK_ID]} -parser:protocol ~/analysis_scripts/match_cm_designs_analysis.xml -parser:view -in:file::extra_res_fa '+params+' -run::preserve_header -enzdes::detect_design_interface -enzdes::cut1 4.0 -enzdes::cut2 6.0 -enzdes::cut3 8.0 -enzdes::cut4 10.0 -enzdes::bb_min_allowed_dev 0.05 -score::weights beta_genpot -corrections::gen_potential -packing::use_input_sc -packing::extrachi_cutoff 1 -packing::ex1 -packing::ex2 -linmem_ig 10 -no_optH false -no_his_his_pairE -lj_hbond_hdis 1.75 -lj_hbond_OH_donor_dis 2.6 -linmem_ig 10 -nblist_autoupdate true -flip_HNQ -out:path:all '+outputdir+' -out:prefix analyze_ -load_PDB_components False -out:file:scorefile '+scorefilename)
    sf.write('\nqstat -j "$JOB_ID"')
    sf.close()
    count+=1
    os.chdir('..')
'''
~/main/source/bin/rosetta_scripts.default.linuxgccrelease -nstruct 1 -jd2:ntrials 1 -s 9UM_1_Q76W45W38_1_model_55706_72512_43281_74068_1_clean_0001.pdb -parser:protocol ~/analysis_scripts/match_cm_designs_analysis.xml -parser:view -in:file::extra_res_fa /wynton/home/kortemme/cgalvin/nps_cm_cst/UM_1_Q76W45W38_1_model_55706_72512_43281_74068_1_lig_ligH_bcc.params -run::preserve_header -enzdes::detect_design_interface -enzdes::cut1 4.0 -enzdes::cut2 6.0 -enzdes::cut3 8.0 -enzdes::cut4 10.0 -enzdes::bb_min_allowed_dev 0.05 -score::weights beta_genpot -corrections::gen_potential -packing::use_input_sc -packing::extrachi_cutoff 1 -packing::ex1 -packing::ex2 -linmem_ig 10 -no_optH false -no_his_his_pairE -lj_hbond_hdis 1.75 -lj_hbond_OH_donor_dis 2.6 -linmem_ig 10 -nblist_autoupdate true -flip_HNQ -out:path:all anl -out:prefix analyze_ -load_PDB_components False -out:file:scorefile test.sc
'''

###
import os
design_dirs=[]
for i in os.listdir():
    if os.path.isdir(i)==True:
        if i[:5]=='match':
            design_dirs.append(i)

for designs_dir in design_dirs:
    os.chdir(designs_dir)
    sfn=[i for i in os.listdir() if i[-3:]=='.sh']
    print(sfn)
    os.system('chmod ugo+x '+ sfn[0])
    os.system('qsub -cwd -t 1-100 -l mem_free=2G '+ sfn[0])
    os.chdir('..')

#

#now consolidate all score files to one directory
import os
design_dirs=[]
for i in os.listdir():
    if os.path.isdir(i)==True:
        if i[:5]=='match':
            design_dirs.append(i)

scoresdir='nps_match2_cm__cst_designs'
os.mkdir(scoresdir)
cwd=os.getcwd()
scoresdirpath=os.path.join(cwd,scoresdir)
for designs_dir in design_dirs:
    os.chdir(designs_dir)
    os.chdir('anl')
    results=[i for i in os.listdir() if i[-3:]=='.sc']
    result=results[0]
    os.system('cp '+result+' '+scoresdirpath+'/'+result)
    os.chdir('..')
    os.chdir('..')

'''
time python ~/tools/merge_all_sf_sameterms.py nps_match2_cm_cst_designs.sc
'''

'''
ANALYSIS
f3=filter(strc_names,datas,'dsasa','>=',0.8)
print(len(f3))

# f5=filter(f3,datas,'ddg','<=',-16.0)
# print(len(f5))


# f6=filter(f5,datas,'hbtolig','>=',1.0)
# print(len(f6))

# f7=filter(f6,datas,'packstat','>=',0.6)
# print(len(f7))

# f8=filter(f7,datas,'boltz','<=',-0.35)
# print(len(f8))

15
['analyze_66_UM_1_Y59I37K57_1_model_293197_63814_24613_20243_1_clean_0001_0001', 'analyze_92_UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_55_UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_61_UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_82_UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_24_UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_57_UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_81_UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_1_UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_70_UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_46_UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_15_UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_39_UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_35_UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_28_UM_1_W55Q68W33_1_model_63022_27975_72512_25825_1_clean_0001_0001']


topstrc 20th perc
[('analyze_89_UM_1_Y27W49Q42_1_model_43674_23554_25797_66239_1_clean_0001_0001', 12),
('analyze_81_UM_1_Q39F60F33_1_model_338872_66238_25402_23019_1_clean_0001_0001', 11),
('analyze_70_UM_1_Q39F60F33_1_model_338872_66238_25402_23019_1_clean_0001_0001', 11),
('analyze_78_UM_1_W39W35L75_1_model_294131_25986_71926_27159_1_clean_0001_0001', 11),
('analyze_10_UM_1_W44W27W29_1_model_148859_74068_65449_18840_1_clean_0001_0001', 10),
('analyze_56_UM_1_Q39F60F33_1_model_338872_66238_25402_23019_1_clean_0001_0001', 10),
('analyze_6_UM_1_Q39F60F33_1_model_338872_66238_25402_23019_1_clean_0001_0001', 10),
('analyze_65_UM_1_Q39F60F33_1_model_338872_66238_25402_23019_1_clean_0001_0001', 10),
('analyze_28_UM_1_Q39F60F33_1_model_338872_66238_25402_23019_1_clean_0001_0001', 10),
('analyze_66_UM_1_W74R81F37_1_model_58134_71926_17596_13854_1_clean_0001_0001', 10)]

ss=[]

import os
design_dirs=[]
for i in os.listdir():
    if os.path.isdir(i)==True:
        if i[:5]=='match':
            design_dirs.append(i)


ss2=[]
for i in ss:
    p=i.split('_')[1:-1]
    n='_'.join(p)
    print(n)
    ss2.append(n+'.pdb')

resultdir=os.path.join(os.getcwd(),'select_strc')
os.mkdir(resultdir)
motif_path='/wynton/home/kortemme/cgalvin/bsff_motifs/nps_2_motif_pdbs/clean_motifs/nps_filtered_motifs/'
for designs_dir in design_dirs:
    os.chdir(designs_dir)
    pdbs=[i for i in os.listdir() if i[-3:]=='pdb']
    for des in ss2:
        if des in pdbs:
            motifres=des.split('_')[7:10]
            motif=motif_path+'_'.join(motifres)+'.pdb'
            os.system('cp '+motif+' '+resultdir+'/'+motif.split('/')[-1])
            os.system('cp '+des+' '+resultdir+'/'+des)
    os.chdir('..')

                cst nps
f=filter(strc_names,datas,'dsasa','>=',0.8)
print(len(f))
846

f2=filter(f,datas,'packstat','>=',0.55)
print(len(f2))
556

f3=filter(f2,datas,'hbtolig','>=',1)
print(len(f3))
223

f4=filter(f3,datas,'boltz','<=',-0.33)
print(len(f4))
f5=filter(f4,datas,'ddg','<=',-15.0)
print(len(f5))
62
61

f6=filter(f5,datas,'total_score','<=',-240)
print(len(f6))
29

f6
ss=['analyze_68UM_1_F29F46Q53_1_model_206131_25403_62132_72512_1_clean_0001_0001', 'analyze_20UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_90UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_27UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_57UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_9UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_19UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_38UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_6UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_1UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_26UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_59UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_36UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_48UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_94UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_8UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_78UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_56UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_35UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_63UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_13UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_67UM_1_W29W40Q1_1_model_327168_25797_18821_63040_1_clean_0001_0001', 'analyze_41UM_1_Y27W49Q42_1_model_43674_23554_25797_66239_1_clean_0001_0001', 'analyze_63UM_1_Y27W49Q42_1_model_43674_23554_25797_66239_1_clean_0001_0001', 'analyze_7UM_1_Y27W49Q42_1_model_43674_23554_25797_66239_1_clean_0001_0001', 'analyze_98UM_1_Y60I37K58_1_model_262691_63814_24613_20243_1_clean_0001_0001', 'analyze_85UM_1_F74F45Q33_1_model_183958_25403_62132_72512_1_clean_0001_0001', 'analyze_89UM_1_F74F45Q33_1_model_183958_25403_62132_72512_1_clean_0001_0001', 'analyze_79UM_1_W55Q68W33_1_model_63022_27975_72512_25825_1_clean_0001_0001']

import os
design_dirs=[]
for i in os.listdir():
    if os.path.isdir(i)==True:
        if i[:5]=='match':
            design_dirs.append(i)

ssdir='npscmcst_select'
os.mkdir(ssdir)
ssdirpath=os.path.join(os.getcwd(),ssdir)
count=0
for dir in design_dirs:
    print(count)
    for i in ss:
        if i+'.pdb' in os.listdir(os.path.join(dir,'anl')):
            count+=1
            os.system('cp '+os.path.join(dir,'anl',i+'.pdb')+' '+os.path.join(ssdirpath,i+'.pdb'))



ANALYZING SELECT DESIGNS FOR MOTIF,SCAFF DIVERSITY
import os
l=os.listdir()
umatches=[]
for des in l:
    s=des.split('_')
    scaffmatch=s[3]
    scaff='_'.join(s[5:7])
    motif='_'.join(s[7:10])
    umatches.append((scaffmatch,scaff,motif))

unique=set(umatches)
In [7]: len(unique)
Out[7]: 6
{('F29F46Q53', 'model_206131', '25403_62132_72512'),
 ('F74F45Q33', 'model_183958', '25403_62132_72512'),
 ('W29W40Q1', 'model_327168', '25797_18821_63040'),
 ('W55Q68W33', 'model_63022', '27975_72512_25825'),
 ('Y27W49Q42', 'model_43674', '23554_25797_66239'),
 ('Y60I37K58', 'model_262691', '63814_24613_20243')}


d={}
for e in unique:
    c=0
    for i in umatches:
        if i==e:
            c+=1
    d[e[0]]=c

{'F74F45Q33': 2,
 'Y60I37K58': 1,
 'W55Q68W33': 1,
 'Y27W49Q42': 3,
 'F29F46Q53': 1,
 'W29W40Q1': 21}


'''




################################################
################################################
################################################
######################OLD#######################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
######################OLD#######################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################






import os
for f in os.listdir():
    if f[-11:-6]=='match':
        cstname=f[:-4]+'.cst'
        os.system('time ipython ~/tools/motif_to_cst.py '+f+' '+paramspath+' '+cstname)

#now to get match pdbs
#in motifs directory that now has matches directory
#make sure no other directories in here
import os
matchpaths=[]
for i in os.listdir():
    if os.path.isdir(i)==True:#directories with motif matche
        os.chdir(i)
        for x in os.listdir():
            if os.path.isdir(x)==True: #directories with scaffolds matches for motif
                os.chdir(x)
                for p in os.listdir():
                    if p[-6:]=='pdb.gz':
                        b=str(os.getcwd())+'/'+p
                        matchpaths.append(b)
                os.chdir('..')
        os.chdir('..')
os.mkdir('matches')
for matchpath in matchpaths:
    fn1=matchpath.split('/')[-3]
    fn2=matchpath.split('/')[-2]
    fn3=matchpath.split('/')[-1][:7]
    if fn3=='scaffol':
        fn=fn1+'_'+fn2+'_'+fn3+'d.pdb.gz'
    else:
        fn=fn1+'_'+fn2+'_'+fn3+'.pdb.gz'
    os.system('cp '+matchpath+' matches/'+fn)

'''
now in matches directory:
gunzip *.gz

matcher needs:
motif cst
scaffold pdb
ligand params
pos file

pos file is just resnumbers to check with 2 spaces after them

motif to cst usage:
#time ipython motif_to_cst.py motif.pdb ligand.params cstfile.cst

'''
#getting motif cst files in matches dir
#I might consider thinking more about constraints in these, make them as strict as possible
paramspath='/wynton/home/kortemme/cgalvin/dexamethasone/dex/Inputs/Rosetta_Inputs/dex.params'
import os
for f in os.listdir():
    if f[-11:-6]=='match':
        cstname=f[:-4]+'.cst'
        os.system('time ipython ~/tools/motif_to_cst.py '+f+' '+paramspath+' '+cstname)

'''
to check a cst file and see if it makes pdb right
~/desktop/Rosetta/main/source/bin/CstfileToTheozymePDB.macosclangrelease -extra_res_fa <.params file of ligand> -match:geometric_constraint_file <cstfile> -database <rosetta database>

scp cgalvin@log2.wynton.ucsf.edu:dexamethasone/motifs/matches/clean_2_17829_21671_match_scaffold_model_396680_match_0.cst test.cst
scp cgalvin@log2.wynton.ucsf.edu:dexamethasone/dex/Inputs/Rosetta_Inputs/dex.params dext.params
~/desktop/Rosetta/main/source/bin/CstfileToTheozymePDB.macosclangrelease -extra_res_fa dex.params -match:geometric_constraint_file test.cst
'''
#get scaffold pos files in matches dir
import os
from pyrosetta import *
init('-ignore_unrecognized_res')
for f in os.listdir():
    if f[-12:]=='scaffold.pdb':
        print(f)
        posname=f[:-4]+'.pos'
        p=pose_from_pdb(f)
        nres=p.total_residue()
        of=open(posname,'w')
        for i in range(nres):
            of.write(str(i+1)+'  ')
        of.close()


'''
scp cgalvin@log2.wynton.ucsf.edu:dexamethasone/motifs/matches/clean_2_17829_21671_match_scaffold_model_396680_match_0.cst ~/desktop/test.cst
scp cgalvin@log2.wynton.ucsf.edu:dexamethasone/motifs/matches/clean_2_17829_21671_match_scaffold_model_396680_scaffold.pos ~/desktop/test.pos
scp cgalvin@log2.wynton.ucsf.edu:dexamethasone/motifs/matches/clean_2_17829_21671_match_scaffold_model_396680_scaffold.pdb ~/desktop/test.pdb

matcher_cmd = ['~/desktop/Rosetta/main/source/bin/match.macosclangrelease',
               '-s', 'test.pdb',
               '-extra_res_fa', 'dex.params',
               '-match:geometric_constraint_file', 'test.cst',
               '-match::scaffold_active_site_residues', 'test.pos',
               '-match:output_format', 'PDB',
               '-match:match_grouper', 'SameSequenceGrouper',
               '-match:consolidate_matches',
               '-match:output_matches_per_group', '1',
               '-use_input_sc',
               '-in:ignore_unrecognized_res',
               '-ex1', '-ex2',
               '-enumerate_ligand_rotamers', 'false',
               '-match::lig_name', 'dex',]

cmd=''
cmd=' '.join(matcher_cmd)

clean_2_17829_21671_match_scaffold_model_343425_match_0.cst
'''
#now to run the standard matcher
paramspath='/wynton/home/kortemme/cgalvin/dexamethasone/dex/Inputs/Rosetta_Inputs/dex.params'
import os
for f in os.listdir():
    if f[-4:]=='.cst':
        scaffoldname=pdbname=f[:-11]+'scaffold.pdb'
        posname=scaffoldname[:-4]+'.pos'
        matcher_cmd = ['~/main/source/bin/match.linuxgccrelease',
                       '-s', scaffoldname,
                       '-extra_res_fa', paramspath,
                       '-match:geometric_constraint_file', f,
                       '-match::scaffold_active_site_residues', posname,
                       '-match:output_format', 'PDB',
                       '-match:match_grouper', 'SameSequenceGrouper',
                       '-match:consolidate_matches',
                       '-match:output_matches_per_group', '1',
                       '-use_input_sc',
                       # '-in:ignore_unrecognized_res',
                       '-ex1', '-ex2','-extrachi_cutoff 0',
                       '-enumerate_ligand_rotamers', 'false',
                       '-match::lig_name', 'dex',] ############SET LIGNAME
        cmd=' '.join(matcher_cmd)
        os.system(cmd)

'''
ALRIGHT IT SEEMS TO HAVE WORKED BUT DIDNT FIND ANY MATCHES
DEFINITELY WANNA KNOW THAT ITS WORKING RIGHT AND THIS IS TO BE EXPECTED
NOT THAT IM BLOWING IT SOMEHOW
'''



#okay, those are just the og scaffold and motif pdbs, not spliced together
#i guess what i actually want to do now is run standard match on these matched
#motifs and scaffs specifically so I can get the match pdbs
#so i wanna put motif and scaffold pdbs in their own folder for conveience
'''
from dex stat motif test
'motifs' dir i submitted for matching contains 3 dir now
clean_2_17829_13638/
clean_1_8292_3566/
clean_2_17829_21671/
are these the 3 motifs which had matches?
    yes it appears so
    and then each directory contains other directories named after scaffold of source match
        each of these contains match info json and 'scaffold.pdb.gz'

scaffold: 2lv8_scaffolds/46/model_221420.pdb.gz
motif: /wynton/home/kortemme/cgalvin/dexamethasone/motifs/clean_1_8292_3566.pdb
num_fast_matches: 1
fast_matched_positions: [[43, 47, 74]]


scaffold: 2lv8_scaffolds/8/model_361984.pdb.gz
motif: /wynton/home/kortemme/cgalvin/dexamethasone/motifs/clean_2_17829_21671.pdb
num_fast_matches: 1
fast_matched_positions: [[33, 69, 38]]


scaffold: 2lv8_scaffolds/42/model_48953.pdb.gz
motif: /wynton/home/kortemme/cgalvin/dexamethasone/motifs/clean_2_17829_21671.pdb
num_fast_matches: 1
fast_matched_positions: [[32, 67, 37]]


scaffold: 2lv8_scaffolds/37/model_169262.pdb.gz
motif: /wynton/home/kortemme/cgalvin/dexamethasone/motifs/clean_2_17829_13638.pdb
num_fast_matches: 1
fast_matched_positions: [[80, 96, 52]]


scaffold: 2lv8_scaffolds/36/model_361428.pdb.gz
motif: /wynton/home/kortemme/cgalvin/dexamethasone/motifs/clean_2_17829_21671.pdb
num_fast_matches: 1
fast_matched_positions: [[33, 69, 38]]


scaffold: 2lv8_scaffolds/12/model_169602.pdb.gz
motif: /wynton/home/kortemme/cgalvin/dexamethasone/motifs/clean_2_17829_13638.pdb
num_fast_matches: 1
fast_matched_positions: [[83, 99, 55]]


scaffold: 2lv8_scaffolds/40/model_60730.pdb.gz
motif: /wynton/home/kortemme/cgalvin/dexamethasone/motifs/clean_2_17829_21671.pdb
num_fast_matches: 1
fast_matched_positions: [[34, 66, 39]]


scaffold: 2lv8_scaffolds/21/model_343425.pdb.gz
motif: /wynton/home/kortemme/cgalvin/dexamethasone/motifs/clean_2_17829_21671.pdb
num_fast_matches: 1
fast_matched_positions: [[36, 78, 41]]


scaffold: 2lv8_scaffolds/7/model_396680.pdb.gz
motif: /wynton/home/kortemme/cgalvin/dexamethasone/motifs/clean_2_17829_21671.pdb
num_fast_matches: 1
fast_matched_positions: [[35, 73, 40]]

'''
