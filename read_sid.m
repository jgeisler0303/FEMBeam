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

function sid= read_sid(fn)
fid= fopen(fn, 'r');
il= 0;

while ~feof(fid)
    line= fgetl(fid);
    if ~isempty(line)
        line= strtrim(line);
    end
    il= il+1;

    if strcmp(line, 'part')
        [sid, il]= read_part(sid, fid, il);
    else
        if ~isempty(line)
            sid.comment{end+1}= line;
        end
    end
end
fclose(fid);


function [sid, il]= read_part(sid, fid, il)
while ~feof(fid)
    line= strtrim(fgetl(fid));
    il= il+1;

    if strncmp(line, 'new modal', 9)
        [sid, il]= read_modal(sid, fid, il);
    else
        switch line
            case 'end part'
                return
            otherwise
                error('no element keyword of class modal found in line %d. found "%s" instead', il, line);
        end
    end
end


function [sid, il]= read_modal(sid, fid, il)
while ~feof(fid)
    line= strtrim(fgetl(fid));
    il= il+1;

    switch line
        case 'end modal'
            return
        case 'refmod'
            [sid.refmod, il]= read_refmod(fid, il);
        case 'frame'
            [sid, il]= read_frame(sid, fid, il);
        case 'mdCM'
            [sid.md, il]= read_taylor(fid, il, 'mdCM');
        case 'J'
            [sid.I, il]= read_taylor(fid, il, 'J');
        case 'Ct'
            [sid.Ct, il]= read_taylor(fid, il, 'Ct');
        case 'Cr'
            [sid.Cr, il]= read_taylor(fid, il, 'Cr');
        case 'Me'
            [sid.Me, il]= read_taylor(fid, il, 'Me');
        case 'Gr'
            [sid.Gr, il]= read_taylor(fid, il, 'Gr');
        case 'Ge'
            [sid.Ge, il]= read_taylor(fid, il, 'Ge');
        case 'Oe'
            [sid.Oe, il]= read_taylor(fid, il, 'Oe');
        case 'ksigma'
            [sid.ksigma, il]= read_taylor(fid, il, 'ksigma');
        case 'Ke'
            [sid.Ke, il]= read_taylor(fid, il, 'Ke');
        case 'De'
            [sid.De, il]= read_taylor(fid, il, 'De');
        otherwise
            error('no element keyword of class modal found in line %d. found "%s" instead', il, line);
    end
end


function i= read_integer(line, il)
i= inf;
pos= strfind(line, '=');
if ~isempty(pos)
    i= str2double(strtrim(line(pos(end)+1:end)));
end
if ~mod(i, 1)==0
    error('no integer found in line %d. found "%s" instead', il, line);
end        

function n= read_number(line, il)
n= nan;
pos= strfind(line, '=');
if ~isempty(pos)
    n= str2double(strrep(strtrim(line(pos(end)+1:end)), 'D', 'e'));
end
if isnan(n)
    error('no number found in line %d. found "%s" instead', il, line);
end        


function m= read_matrix(m, line, il, ni)
pos1= strfind(line, '(');
pos2= strfind(line, ')');
if isempty(pos1) || isempty(pos2) || pos1(end)>=pos2(1)
    error('no matrix index found in line %d. found "%s" instead', il, line)
end
idx= str2num(line(pos1(end)+1:pos2(1)-1));
if ~all(mod(idx, 1)==0) || length(idx)~=ni
    error('no proper matrix index found in line %d. found "%s" instead', il, line)
end
if ni==2
    m(idx(1), idx(2))= read_number(line, il);
else
    m(idx(1), idx(2), idx(3))= read_number(line, il);
end    


function c= read_cell(c, line, il)
pos1= strfind(line, '(');
pos2= strfind(line, ')');
if isempty(pos1) || isempty(pos2) || pos1(end)>=pos2(1)
    error('no cell index found in line %d. found "%s" instead', il, line)
end
idx= str2num(line(pos1(end)+1:pos2(1)-1));
if ~all(mod(idx, 1)==0) || length(idx)~=1
    error('no proper cell index found in line %d. found "%s" instead', il, line)
end
pos= strfind(line, '=');
if ~isempty(pos)
    c{idx(1)}= strtrim(line(pos(end)+1:end));
else
    c{idx(1)}= 'noname';
end


function [t, il]= read_taylor(fid, il, name)
while ~feof(fid)
    line= strtrim(fgetl(fid));
    il= il+1;

    if strncmpi(lower(line), 'm0', 2)
        if ~isfield(t, 'M0')
            if ~isfield(t, 'order') || ~isfield(t, 'nrow') || ~isfield(t, 'ncol') || ~isfield(t, 'nq') || ~isfield(t, 'nqn') || ~isfield(t, 'structure')
                error('found m0 entry but not all of order, nrow, ncol, np, npn, and structure have been defined in line %d. found "%s" instead', il, line);
            end
            t.M0= zeros(t.nrow, t.ncol);
        end
        t.M0= read_matrix(t.M0, line, il, 2);
    elseif strncmpi(lower(line), 'm1', 2)
        if ~isfield(t, 'M1')
            if ~isfield(t, 'order') || ~isfield(t, 'nrow') || ~isfield(t, 'ncol') || ~isfield(t, 'nq') || ~isfield(t, 'nqn') || ~isfield(t, 'structure')
                error('found m1 entry but not all of order, nrow, ncol, np, npn, and structure have been defined in line %d. found "%s" instead', il, line);
            end
            t.M1= zeros(t.nrow, t.ncol, t.nq);
        end
        t.M1= read_matrix(t.M1, line, il, 3);
    elseif strncmpi(line, 'mn', 2)
        if ~isfield(t, 'mn')
            if ~isfield(t, 'order') || ~isfield(t, 'nrow') || ~isfield(t, 'ncol') || ~isfield(t, 'nq') || ~isfield(t, 'nqn') || ~isfield(t, 'structure')
                error('found mn entry but not all of order, nrow, ncol, np, npn, and structure have been defined in line %d. found "%s" instead', il, line);
            end
            t.mn= zeros(t.nrow, t.ncol, t.nqn);
        end
        t.mn= read_matrix(t.mn, line, il, 3);
    elseif strncmpi(line, 'order', 5)
        t.order= read_integer(line, il);
    elseif strncmpi(line, 'nrow', 4)
        t.nrow= read_integer(line, il);
    elseif strncmpi(line, 'ncol', 4)
        t.ncol= read_integer(line, il);
    elseif strncmpi(line, 'nqn', 3)
        t.nqn= read_integer(line, il);
    elseif strncmpi(line, 'nq', 2)
        t.nq= read_integer(line, il);
    elseif strncmpi(line, 'structure', 6)
        t.structure= read_integer(line, il);
    elseif strcmp(lower(line), ['end ' lower(name)])
        if ~isfield(t, 'order') || ~isfield(t, 'nrow') || ~isfield(t, 'ncol') || ~isfield(t, 'nq') || ~isfield(t, 'structure')
            error('taylor element ends but order, nrow, ncol, structure, and nq have not been defined in line %d', il);
        elseif t.order>1 && ~isfield(t, 'nqn') 
            error('taylor element ends but nqn has not been defined in line %d', il);
        else
            if ~isfield(t, 'M0')
                t.M0= zeros(t.nrow, t.ncol);
            end
            if t.order>0 && ~isfield(t, 'M1')
                t.M1= zeros(t.nrow, t.ncol, t.nq);
            end
            if t.order>1 && ~isfield(t, 'mn')
                t.mn= zeros(t.nrow, t.ncol, t.nqn);
            end            
        end
        return
    else
        error('no element keyword of class taylor found in line %d. found "%s" instead', il, line);
    end
end


function [node, il]= read_node(fid, il, name)
node.node= name;
while ~feof(fid)
    line= strtrim(fgetl(fid));
    il= il+1;

    if strncmp(line, 'rframe', 6)
        pos= strfind(line, '=');
        if ~isempty(pos)
            node.rframe= strtrim(line(pos(end)+1:end));
        else
            node.rframe= 'unknown ref';
        end        
    else
        switch lower(line)
            case 'end node'
                return
            case 'origin'
                [node.origin, il]= read_taylor(fid, il, 'origin');
            case 'ap'
                [node.AP, il]= read_taylor(fid, il, 'AP');
            case 'phi'
                [node.Phi, il]= read_taylor(fid, il, 'Phi');
            case 'psi'
                [node.Psi, il]= read_taylor(fid, il, 'Psi');
            case 'sigma'
                [node.sigma, il]= read_taylor(fid, il, 'sigma');
            otherwise
                error('no element keyword of class node found in line %d. found "%s" instead', il, line);
        end
    end
end


function [sid, il]= read_frame(sid, fid, il)
while ~feof(fid)
    line= strtrim(fgetl(fid));
    il= il+1;

    if strcmp(line, 'end frame')
        return
    end
    if strncmp(line, 'new node', 8)
        pos= strfind(line, '=');
        if ~isempty(pos)
            name= strtrim(line(pos(end)+1:end));
        else
            name= 'noname';
        end    
        [sid.frame(end+1), il]= read_node(fid, il, name);
    else
        error('no element keyword of class frame found in line %d. found "%s" instead', il, line);
    end
end

function [refmod, il]= read_refmod(fid, il)
refmod= struct();
while ~feof(fid)
    line= strtrim(fgetl(fid));
    il= il+1;

    if strncmpi(line, 'mass', 4)
        refmod.mass= read_number(line, il);
    elseif strncmpi(line, 'nelastq', 7)
        refmod.nelastq= read_integer(line, il);
        refmod.ielastq= cell(refmod.nelastq, 1);
    elseif strncmpi(line, 'ielastq', 7)
        refmod.ielastq= read_cell(refmod.ielastq, line, il);
    elseif strcmp(line, 'end refmod')
        return
    else
        error('no element keyword of class refmod refmod in line %d. found "%s" instead', il, line);
    end
end



