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
 

function write_sid(sid, fn, nodes)
[pn, fn, ex]= fileparts(fn);
if isempty(ex)
    ex= '.SID_FEM';
end
fn= fullfile(pn, [fn ex]);

fid= fopen(fn, 'w');

com= sid.comment;
if ~iscell(com)
    com= {com};
end
for i= 1:length(com)
    fprintf(fid, '%s\n', com{i});
end
fprintf(fid, 'part\n');
fprintf(fid, '    new modal\n'); % what to write when body is not modal?
write_refmod(sid.refmod, fid, '        ');

fprintf(fid, '        frame\n'); % what to write when body is not modal?
for i= 1:length(sid.frame)
    if exist('nodes', 'var') && ~ismember(i, nodes), continue; end
        
    write_frame(sid.frame(i), fid, '                ');
end    
fprintf(fid, '        end frame\n');    

write_taylor('mdCM', sid.mdCM, fid, '        ');
write_taylor('J', sid.J, fid, '        ');
write_taylor('Ct', sid.Ct, fid, '        ');
write_taylor('Cr', sid.Cr, fid, '        ');
write_taylor('Me', sid.Me, fid, '        ');
write_taylor('Gr', sid.Gr, fid, '        ');
write_taylor('Ge', sid.Ge, fid, '        ');
write_taylor('Oe', sid.Oe, fid, '        ');
write_taylor('ksigma', sid.ksigma, fid, '        ');
write_taylor('Ke', sid.Ke, fid, '        ');
write_taylor('De', sid.De, fid, '        ');

fprintf(fid, '    end modal\n');
fprintf(fid, 'end part\n');

fclose(fid);

function write_refmod(refmod, fid, indent)
fprintf(fid, '%srefmod\n', indent);
fprintf(fid, '%s    mass           = %f\n', indent, refmod.mass);
fprintf(fid, '%s    nelastq        = %d\n', indent, refmod.nelastq);
for i= 1:length(refmod.ielastq)
    fprintf(fid, '%s    ielastq (%4d) = %s\n', indent, i, refmod.ielastq{i});
end
fprintf(fid, '%send refmod\n', indent);


function write_frame(f, fid, indent)
fprintf(fid, '%snew node   = %s\n', indent, f.node);

fprintf(fid, '%s    rframe = %s\n', indent, f.rframe);
write_taylor('origin', f.origin, fid, [indent '    ']);
write_taylor('phi', f.Phi, fid, [indent '    ']);
write_taylor('psi', f.Psi, fid, [indent '    ']);
write_taylor('AP', f.AP, fid, [indent '    ']);
if isfield(f, 'sigma')
    write_taylor('sigma', f.sigma, fid, [indent '    ']);
end

fprintf(fid, '%send node\n', indent);



function write_taylor(name, t, fid, indent)
fprintf(fid, '%s%s\n', indent, name);
fprintf(fid, '%s    order        = %d\n', indent, t.order);
fprintf(fid, '%s    nrow         = %d\n', indent, t.nrow);
fprintf(fid, '%s    ncol         = %d\n', indent, t.ncol);
fprintf(fid, '%s    nq           = %d\n', indent, t.nq);
fprintf(fid, '%s    nqn          = %d\n', indent, t.nqn);
fprintf(fid, '%s    structure    = %d\n', indent, t.structure);
structure= t.structure;
if isempty(t.M0), structure= 0; end
switch structure
    case 0
    case 1
        for i= 1:t.nrow
            fprintf(fid, '%s    m0(%2d,%2d)    = %f\n', indent, i, i, t.M0(i, i));
        end
    case 2
        for i= 1:t.nrow
            for j= 1:i
                fprintf(fid, '%s    m0(%2d,%2d)    = %f\n', indent, i, j, t.M0(i, j));
            end
        end
    case 3
        for i= 1:t.nrow
            for j= 1:t.ncol
                fprintf(fid, '%s    m0(%2d,%2d)    = %f\n', indent, i, j, t.M0(i, j));
            end
        end
    case 4
end

if t.order>0
    for i= 1:t.nrow
        for j= 1:t.ncol
            for k= 1:t.nq
                fprintf(fid, '%s    m1(%2d,%2d,%2d) = %f\n', indent, i, j, k, t.M1(i, j, k));
            end
        end
    end
end

fprintf(fid, '%send %s\n', indent, name);
