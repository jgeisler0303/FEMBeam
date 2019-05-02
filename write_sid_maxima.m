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
    fprintf(fid, '%s@ielastq[%d]: "%s";\n', body_name, i, sid.refmod.ielastq{i});
end
fprintf(fid, '\n');

% TODO: write frame data

write_taylor(fid, body_name, 'md', sid.md, false, threshold);
write_taylor(fid, body_name, 'I', sid.I, false, threshold);
write_taylor(fid, body_name, 'Ct', sid.Ct, false, threshold);
write_taylor(fid, body_name, 'Cr', sid.Cr, false, threshold);
write_taylor(fid, body_name, 'Me', sid.Me, false, threshold);
write_taylor(fid, body_name, 'Gr', sid.Gr, true, threshold);
write_taylor(fid, body_name, 'Ge', sid.Ge, true, threshold);
write_taylor(fid, body_name, 'Oe', sid.Oe, false, threshold);
%  write_taylor(fid, body_name, 'ksigma', sid.ksigma, threshold);
write_taylor(fid, body_name, 'Ke', sid.Ke, false, threshold);
write_taylor(fid, body_name, 'De', sid.De, false, threshold);

fclose(fid);

%  function write_frame(f, fid, indent)
%  fprintf(fid, '%snew node   = %s\n', indent, f.node);
%  
%  fprintf(fid, '%s    rframe = %s\n', indent, f.rframe);
%  write_taylor('origin', f.origin, fid, [indent '    ']);
%  write_taylor(fid, 'phi', f.Phi, [indent '    ']);
%  write_taylor(fid, 'psi', f.Psi, [indent '    ']);
%  write_taylor(fid, 'AP', f.AP, [indent '    ']);
%  if isfield(f, 'sigma')
%      write_taylor(fid, 'sigma', f.sigma, [indent '    ']);
%  end
%  
%  fprintf(fid, '%send node\n', indent);


function write_taylor(fid, name, var_name, t, with_sub_idx, threshold)
structure= t.structure;
if isempty(t.M0), structure= 0; end
    
if with_sub_idx
    nelem= 3;
else
    nelem= 0;
end
switch structure
    case 0
    case 1
        for i= 1:t.nrow
            write_t_element(fid, sprintf('%s@%s', name, var_name), t.M0, [i, i], nelem, threshold);
        end
    case {2, 3}
        for i= 1:t.nrow
            for j= 1:t.ncol
                write_t_element(fid, sprintf('%s@%s', name, var_name), t.M0, [i, j], nelem, threshold);
            end
        end
    case 4
end

if t.order>0
    for k= 1:t.nq
        for i= 1:t.nrow
            for j= 1:t.ncol
                write_t_element(fid, sprintf('%s@%s', name, var_name), t.M1, [i, j, k], nelem, threshold);
            end
        end
    end
end
fprintf(fid, '\n');

function write_t_element(fid, var_name, data, idx, nelem, threshold)
v= subsref(data, struct('type', '()', 'subs', {num2cell(idx)}));
if abs(v)<threshold, return; end
    
if length(idx)==2
    write_m_element(fid, sprintf('%s@M0', var_name), v, idx, nelem);
else
    write_m_element(fid, sprintf('%s@M1[%d]', var_name, idx(3)), v, idx(1:2), nelem);
end

function write_m_element(fid, var_name, v, idx, nelem)
if nelem==0
    fprintf(fid, '%s[%d, %d]: %g\n', var_name, idx(1), idx(2), v);
else
    fprintf(fid, '%s[%d][%d, %d]: %g\n', var_name, fix((idx(2)-1)/nelem)+1, idx(1), rem((idx(2)-1), nelem)+1, v);
end

