function model = comsolMSR(pos)
%
% SNS_MSR_Chris_S_EDIT.m
%
% Model exported on Nov 7 2017, 18:01 by COMSOL 5.3.0.260.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('/data1/cmswank/SpinSimSwank2/data');

model.label('SNS_MSR_Chris_S_EDIT.mph');

model.comments(['Untitled\n\n']);

model.param.set('Cur', '2.2169 [A]');
model.param.set('Yhalf', '2 [m]');
model.param.set('Zhalf', '3 [m]');
model.param.set('coilGap', '20 [cm]');
model.param.set('ShieldThick', '5 [mm]');
model.param.set('interthick', '5[mm]');
model.param.set('shieldcoilgap', '10 [cm]');

model.component.create('comp1', false);

model.component('comp1').geom.create('geom1', 3);



model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').label('Magnetic_Setup');
model.component('comp1').geom('geom1').geomRep('comsol');
model.component('comp1').geom('geom1').repairTolType('relative');
model.component('comp1').geom('geom1').selection.create('csel1', 'CumulativeSelection');
model.component('comp1').geom('geom1').selection('csel1').label('Coils');
model.component('comp1').geom('geom1').selection.create('csel2', 'CumulativeSelection');
model.component('comp1').geom('geom1').selection('csel2').label('Shield');
model.component('comp1').geom('geom1').selection.create('csel3', 'CumulativeSelection');
model.component('comp1').geom('geom1').selection('csel3').label('AirMat');
model.component('comp1').geom('geom1').selection.create('csel4', 'CumulativeSelection');
model.component('comp1').geom('geom1').selection('csel4').label('CoilsRotate');
model.component('comp1').geom('geom1').selection.create('csel5', 'CumulativeSelection');
model.component('comp1').geom('geom1').selection('csel5').label('Sim Bound');
model.component('comp1').geom('geom1').create('blk3', 'Block');
model.component('comp1').geom('geom1').feature('blk3').set('contributeto', 'csel2');
model.component('comp1').geom('geom1').feature('blk3').set('base', 'center');
model.component('comp1').geom('geom1').feature('blk3').set('size', {'2*(Yhalf + shieldcoilgap)' '2*(Yhalf + shieldcoilgap)' '2*(Zhalf + shieldcoilgap)'});
model.component('comp1').geom('geom1').create('blk4', 'Block');
model.component('comp1').geom('geom1').feature('blk4').set('contributeto', 'csel5');
model.component('comp1').geom('geom1').feature('blk4').set('base', 'center');
model.component('comp1').geom('geom1').feature('blk4').set('size', {'2*(Yhalf + shieldcoilgap)*1.2' '2*(Yhalf + shieldcoilgap)*1.2' '2*(Zhalf + shieldcoilgap)*1.2'});

%%% COIL GEOMETRY

model.component('comp1').geom('geom1').create('ic1', 'InterpolationCurve');
model.component('comp1').geom('geom1').feature('ic1').set('contributeto', 'csel1');
model.component('comp1').geom('geom1').feature('ic1').set('table', [0 0 0; 1 1 1]);
model.component('comp1').geom('geom1').create('ic2', 'InterpolationCurve');
model.component('comp1').geom('geom1').feature('ic2').set('contributeto', 'csel1');
model.component('comp1').geom('geom1').feature('ic2').set('table', [1 1 1; 1 1 0]);
model.component('comp1').geom('geom1').create('ic3', 'InterpolationCurve');
model.component('comp1').geom('geom1').feature('ic3').set('contributeto', 'csel1');
model.component('comp1').geom('geom1').feature('ic3').set('table', [1 1 0; 0 0 0]);
model.component('comp1').geom('geom1').feature('fin').set('repairtoltype', 'relative');
model.component('comp1').geom('geom1').run;



model.component('comp1').coordSystem.create('sys2', 'Cylindrical');

model.component('comp1').physics.create('mf', 'InductionCurrents', 'geom1');
model.component('comp1').physics('mf').create('edc1', 'EdgeCurrent', 1);
model.component('comp1').physics('mf').feature('edc1').selection.named('geom1_csel1_edg');
model.component('comp1').physics('mf').create('ms2', 'MagneticShielding', 2);
model.component('comp1').physics('mf').feature('ms2').selection.named('geom1_csel2_bnd');

model.component('comp1').mesh('mesh1').autoMeshSize(4);




model.component('comp1').coordSystem('sys1').set('mastercoordsystcomp', 'manual');
model.component('comp1').coordSystem('sys2').set('name', 'cylsys');
model.component('comp1').coordSystem('sys2').set('frametype', 'geometry');

model.component('comp1').physics('mf').prop('MeshControl').set('EnableMeshControl', true);
model.component('comp1').physics('mf').feature('edc1').set('Ie', 1);
model.component('comp1').physics('mf').feature('ms2').set('murbnd', 60000);
model.component('comp1').physics('mf').feature('al1').set('sigma_mat', 'userdef');
model.component('comp1').physics('mf').feature('al1').set('epsilonr_mat', 'userdef');
model.component('comp1').physics('mf').feature('al1').set('mur_mat', 'userdef');
model.component('comp1').physics('mf').feature('ms2').set('murbnd_mat', 'userdef');

model.study.create('std1');
model.study('std1').create('stat', 'Stationary');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').create('i1', 'Iterative');
model.sol('sol1').feature('s1').feature('i1').create('mg1', 'Multigrid');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').create('so1', 'SOR');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').create('so1', 'SOR');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').create('ams1', 'AMS');
model.sol('sol1').feature('s1').feature.remove('fcDef');



model.sol('sol1').attach('std1');
model.sol('sol1').feature('s1').feature('i1').set('linsolver', 'fgmres');
model.sol('sol1').runAll;


out = model;
