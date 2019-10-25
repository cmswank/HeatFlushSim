function out = HeatFlushTemplate()
%
% HeatFlushTemplate.m
%
% Model exported on Oct 5 2017, 11:57 by COMSOL 5.3.0.260.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('/data1/cmswank/SpinSimSwank2/data');

model.label('HeatFlush001.mph');

model.comments(['Untitled\n\n']);

model.baseSystem('cgs');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 3);

model.result.table.create('evl3', 'Table');
model.result.table.create('tbl1', 'Table');
model.result.table.create('tbl2', 'Table');

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').create('cyl1', 'Cylinder');
model.component('comp1').geom('geom1').feature('cyl1').set('r', 2.03);
model.component('comp1').geom('geom1').feature('cyl1').set('h', 100);
model.component('comp1').geom('geom1').create('blk1', 'Block');
model.component('comp1').geom('geom1').feature('blk1').set('pos', {'0' '0' '100+7.6/2'});
model.component('comp1').geom('geom1').feature('blk1').set('base', 'center');
model.component('comp1').geom('geom1').feature('blk1').set('size', [40 10.2 7.6]);
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('fin');

model.component('comp1').variable.create('var1');
model.component('comp1').variable('var1').set('Temp', '(p>0)*(p/2492.2)^0.25');
model.component('comp1').variable('var1').set('lam', '1.28E-5*((abs(p)/2492.2)^0.25)^(-9)+3.80E-3*((abs(p)/2492.2)^0.25)^(-4)');
model.component('comp1').variable('var1').set('rhon', '4.43E-7/(.4)^4*Temp^4');
model.component('comp1').variable('var1').set('vth', '2.4e4');

model.component('comp1').physics.create('spf', 'LaminarFlow', 'geom1');
model.component('comp1').physics('spf').create('inl1', 'InletBoundary', 2);
model.component('comp1').physics('spf').feature('inl1').selection.set([8]);
model.component('comp1').physics('spf').create('out1', 'OutletBoundary', 2);
model.component('comp1').physics('spf').feature('out1').selection.set([4]);

model.component('comp1').mesh('mesh1').autoMeshSize(4);

model.result.table('evl3').label('Evaluation 3D');
model.result.table('evl3').comments('Interactive 3D values');
model.result.table('tbl1').comments('Surface Integration 1 (4*w*p)');
model.result.table('tbl2').comments('Surface Average 1 (4*w*p)');

model.component('comp1').view('view1').set('renderwireframe', true);
model.component('comp1').view('view1').set('scenelight', false);

model.component('comp1').physics('spf').prop('PhysicalModelProperty').set('Tref', 0.45);
model.component('comp1').physics('spf').prop('PhysicalModelProperty').set('pref', '1E-2');
model.component('comp1').physics('spf').feature('fp1').set('rho', 'rhon');
model.component('comp1').physics('spf').feature('fp1').set('mu', 'rhon*vth*lam/3');
model.component('comp1').physics('spf').feature('init1').set('u_init', [0; 0; 10]);
model.component('comp1').physics('spf').feature('init1').set('p_init', '(.45/2492.2)^0.25');
model.component('comp1').physics('spf').feature('inl1').set('BoundaryCondition', 'Pressure');
model.component('comp1').physics('spf').feature('inl1').set('U0in', '0.010*1E7/4/p/4.1209');
model.component('comp1').physics('spf').feature('inl1').set('p0', '2492.2*(0.455)^4');
model.component('comp1').physics('spf').feature('out1').set('p0', '2492.2*0.45^4');
model.component('comp1').physics('spf').feature('out1').set('NormalFlow', true);

model.frame('mesh1').sorder(1);

model.component('comp1').physics('spf').feature('fp1').set('rho_mat', 'userdef');
model.component('comp1').physics('spf').feature('fp1').set('mu_mat', 'userdef');

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
model.sol('sol1').feature('s1').create('i2', 'Iterative');
model.sol('sol1').feature('s1').feature('i1').create('mg1', 'Multigrid');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').create('sc1', 'SCGS');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').create('sc1', 'SCGS');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol1').feature('s1').feature('i2').create('mg1', 'Multigrid');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('pr').create('sc1', 'SCGS');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('po').create('sc1', 'SCGS');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol1').feature('s1').feature.remove('fcDef');

model.result.dataset.create('surf1', 'Surface');
model.result.dataset.create('surf2', 'Surface');
model.result.dataset.create('surf3', 'Surface');
model.result.dataset.create('surf4', 'Surface');
model.result.dataset('surf1').selection.set([1 2 3 5 6 7 8 10 11 12]);
model.result.dataset('surf2').selection.set([1 2 3 5 6 7 8 10 11 12]);
model.result.dataset('surf3').selection.set([1 2 3 5 6 7 8 10 11 12]);
model.result.dataset('surf4').selection.set([1 2 3 5 6 7 8 10 11 12]);
model.result.numerical.create('int1', 'IntSurface');
model.result.numerical.create('av1', 'AvSurface');
model.result.numerical('int1').selection.set([9]);
model.result.numerical('int1').set('probetag', 'none');
model.result.numerical('av1').selection.set([9]);
model.result.numerical('av1').set('probetag', 'none');
model.result.create('pg1', 'PlotGroup3D');
model.result('pg1').create('slc1', 'Slice');
model.result('pg1').create('iso1', 'Isosurface');
model.result('pg1').create('str1', 'Streamline');
model.result('pg1').feature('str1').selection.set([8]);

model.sol('sol1').attach('std1');
model.sol('sol1').feature('st1').set('splitcomplex', true);
model.sol('sol1').feature('s1').set('nonlin', true);
model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'i1');
model.sol('sol1').feature('s1').feature('fc1').set('initstep', 0.01);
model.sol('sol1').feature('s1').feature('fc1').set('minstep', 1.0E-6);
model.sol('sol1').feature('s1').feature('fc1').set('maxiter', 100);
model.sol('sol1').feature('s1').feature('i1').label('Algebraic Multigrid Solver (spf)');
model.sol('sol1').feature('s1').feature('i1').set('nlinnormuse', true);
model.sol('sol1').feature('s1').feature('i1').set('maxlinit', 200);
model.sol('sol1').feature('s1').feature('i1').set('rhob', 20);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('prefun', 'saamg');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('mgcycle', 'f');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('maxcoarsedof', 30000);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('strconn', 0.02);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('saamgcompwise', true);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('usesmooth', false);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sc1').set('linesweeptype', 'ssor');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sc1').set('iter', 0);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sc1').set('scgsvertexrelax', 0.7);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sc1').set('linesweeptype', 'ssor');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sc1').set('iter', 1);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sc1').set('scgsvertexrelax', 0.7);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('s1').feature('i2').label('Geometric Multigrid Solver (spf)');
model.sol('sol1').feature('s1').feature('i2').set('nlinnormuse', true);
model.sol('sol1').feature('s1').feature('i2').set('maxlinit', 200);
model.sol('sol1').feature('s1').feature('i2').set('rhob', 20);
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('pr').feature('sc1').set('linesweeptype', 'ssor');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('pr').feature('sc1').set('iter', 0);
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('pr').feature('sc1').set('scgsvertexrelax', 0.7);
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('po').feature('sc1').set('linesweeptype', 'ssor');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('po').feature('sc1').set('iter', 1);
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('po').feature('sc1').set('scgsvertexrelax', 0.7);
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').runAll;

model.result.dataset('surf1').label('Exterior Walls');
model.result.dataset('surf2').label('Exterior Walls 1');
model.result.dataset('surf3').label('Exterior Walls 2');
model.result.dataset('surf4').label('Exterior Walls 3');
model.result.numerical('int1').set('table', 'tbl1');
model.result.numerical('int1').set('expr', {'4*w*p'});
model.result.numerical('int1').set('unit', {'g*cm^2/s^3'});
model.result.numerical('int1').set('descr', {''});
model.result.numerical('av1').set('table', 'tbl2');
model.result.numerical('av1').set('expr', {'4*w*p'});
model.result.numerical('av1').set('unit', {'g/s^3'});
model.result.numerical('av1').set('descr', {''});
model.result.numerical('int1').setResult;
model.result.numerical('av1').setResult;
model.result('pg1').set('legendcolor', 'custom');
model.result('pg1').feature('slc1').set('expr', 'w');
model.result('pg1').feature('slc1').set('descr', 'Velocity field, z component');
model.result('pg1').feature('slc1').set('resolution', 'normal');
model.result('pg1').feature('iso1').active(false);
model.result('pg1').feature('iso1').set('expr', 'spf.U*(x>0)');
model.result('pg1').feature('iso1').set('descr', 'spf.U*(x>0)');
model.result('pg1').feature('iso1').set('interactive', true);
model.result('pg1').feature('iso1').set('shift', -14);
model.result('pg1').feature('iso1').set('resolution', 'normal');
model.result('pg1').feature('str1').set('selnumber', 10);
model.result('pg1').feature('str1').set('linetype', 'tube');
model.result('pg1').feature('str1').set('radiusexpr', 'w');
model.result('pg1').feature('str1').set('tuberadiusscale', 0.006499371927629146);
model.result('pg1').feature('str1').set('tuberadiusscaleactive', false);
model.result('pg1').feature('str1').set('resolution', 'normal');

out = model;
