%  Copyright 2019, Jens Geisler
%  
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%  
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 

function write_sid_maxima(sid, fn, body_name, nodes, threshold)
if ~exist('threshold', 'var')
    threshold= eps;
end

[pn, fn, ex]= fileparts(fn);
if isempty(ex)
    ex= '.mac';
end
fn= fullfile(pn, [fn ex]);

fid= fopen(fn, 'w');

com= sid.comment;
if ~iscell(com)
    com= {com};
end
for i= 1:length(com)
    fprintf(fid, '/* %s */\n', com{i});
end

fprintf(fid, '\n%s: emptyElasticMode(%d);\n\n', body_name, sid.refmod.nelastq);

fprintf(fid, '%s@mass: %g;\n\n', body_name, sid.refmod.nelastq);

for i= 1:sid.refmod.nelastq
    fprintf(fid, '%rm@ielastq[%d]: %s;\n', body_name, i, sid.refmod.ielastq{i});
end

% TODO: write frame data

write_taylor(body_name, 'md', sid.md, false, fid, threshold);
write_taylor(body_name, 'I', sid.I, false, fid, threshold);
write_taylor(body_name, 'Ct', sid.Ct, false, fid, threshold);
write_taylor(body_name, 'Cr', sid.Cr, false, fid, threshold);
write_taylor(body_name, 'Me', sid.Me, false, fid, threshold);
write_taylor(body_name, 'Gr', sid.Gr, true, fid, threshold);
write_taylor(body_name, 'Ge', sid.Ge, true, fid, threshold);
write_taylor(body_name, 'Oe', sid.Oe, false, fid, threshold);
%  write_taylor(body_name, 'ksigma', sid.ksigma, fid, threshold);
write_taylor(body_name, 'Ke', sid.Ke, false, fid, threshold);
write_taylor(body_name, 'De', sid.De, false, fid, threshold);

fclose(fid);

%  function write_frame(f, fid, indent)
%  fprintf(fid, '%snew node   = %s\n', indent, f.node);
%  
%  fprintf(fid, '%s    rframe = %s\n', indent, f.rframe);
%  write_taylor('origin', f.origin, fid, [indent '    ']);
%  write_taylor('phi', f.Phi, fid, [indent '    ']);
%  write_taylor('psi', f.Psi, fid, [indent '    ']);
%  write_taylor('AP', f.AP, fid, [indent '    ']);
%  if isfield(f, 'sigma')
%      write_taylor('sigma', f.sigma, fid, [indent '    ']);
%  end
%  
%  fprintf(fid, '%send node\n', indent);


function write_taylor(name, var_name, t, with_sub_idx, fid, threshold)
structure= t.structure;
if isempty(t.M0), structure= 0; end
    
if with_sub_idx
    nelem= t.nq;
else
    nelem= 0;
end
switch structure
    case 0
    case 1
        for i= 1:t.nrow
            write_m_element(sprintf('%s@%s', name, var_name), t.M0, [i, i], nelem, threshold);
        end
    case {2, 3}
        for i= 1:t.nrow
            for j= 1:t.ncol
                write_m_element(sprintf('%s@%s', name, var_name), t.M0, [i, j], nelem, threshold);
            end
        end
    case 4
end

if t.order>0
    for i= 1:t.nrow
        for j= 1:t.ncol
            for k= 1:t.nq
                write_m_element(sprintf('%s@%s', name, var_name), t.M1, [i, j, k], nelem, threshold);
            end
        end
    end
end


function write_m_element(var_name, data, idx, nelem, threshold)
v= subsref(data, struct('type', '()', 'subs', {num2cell(idx)}));
if abs(v)<threshold, return; end
    
if length(idx)==2
    if nelem>0
        fprintf(fid, '%s@M0[%d, %d]: %g\n', var_name, idx(1), idx(2), v);
    else
        fprintf(fid, '%s@M0[%d, %d]: %g\n', var_name, idx(1), idx(2), v);
    end
else
    fprintf(fid, '%s@M1[%d][%d, %d]: %g\n', var_name, idx(3), idx(1), idx(2), v)
end

