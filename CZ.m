classdef CZ %< handle
    % 结构：n=[n1,n2,...,nk]是群的结构 (double array)
    % fhat=CZ_2n.c fhat(g)就是 CZ_2n的系数 既： g =char array
    properties
        n
        c
    end
    methods
        %% 构造函数
        function obj = CZ(n)
            obj.n=n;
            obj.c=containers.Map('KeyType',  'char', 'ValueType', 'any');
        end
        %% Get/Set
        function  x = get(f,t)
            %得到 函数的系数
            if ~ischar(t)
                t=char(t);
            end
            if length(t(:))==length(f.n)
                t=t(:)';
                if isKey(f.c,t)
                    x=f.c(t);
                else
                    x=0;
                end
            end
        end
        
        function  f = set(f,t,x)
            t=double(t);
            t=mod(t,f.n);
            t=char(t);
            f.c(t)=x;
        end
        %% 赋值
        function x = subsref(this,s)
            if strcmp(s.type,'{}')
                x=cell2mat(s.subs);
                x=PointValue(this,x);
                return
            end
            if strcmp(s.subs,'c')&&strcmp(s.type,'.')
                x=this.c;
                return
            end
            if strcmp(s.subs,'n')&&strcmp(s.type,'.')
                x=this.n;
                return
            end
            s=s.subs;
            t=s{1};
            if ~ischar(t)
                t=mod(t,this.n);
                t=char(t);
            end
            if isKey(this.c,t)
                x=get(this,t);
            else
                x=0;
            end
        end
        
        function this = subsasgn(this,s,b)
            s=s.subs;
            t=s{1};
            if ~ischar(t)
                t=char(t);
            end
            set(this,t,b);
        end
        
        function this=setcoeff(this,c)
            this.c=c;
        end
        %% 加法
        function h=plus(f,g)
            if isnumeric(f)&&length(f)==1
                h=g;
                N=length(h.n);t=zeros(1,N);
                set(h,t,get(h,t)+f);
                return
            end
            if isnumeric(g)&&length(g)==1
                h=plus(g,f);
                return
            end
            if ~(all(f.n==g.n))
                error('the addition are not on same group')
            end
            h=CZ(f.n);
            Index=[cell2mat(keys(f.c).');cell2mat(keys(g.c).')];
            Index=unique(Index,'rows');
            for i=1:(size(Index,1))
                t=Index(i,:);
                if get(f,t)+get(g,t)~=0
                    set(h,t,get(f,t)+get(g,t));
                end
            end
        end
        %% Sparsity
        function L=length(this)
            L=length(this.c);
        end
        function L=sparsity(this)
            L=length(this.c);
        end
        %% 减法
        function h = minus(f,g)
            if isnumeric(f)&&length(f)==1
                h=minus(CZ(g.n),g)+f;
                return
            end
            if isnumeric(g)&&length(g)==1
                h=plus(f,-g);
                return
            end
            if f.n~=g.n
                error('the addition are not on same group')
            end
            h=CZ(f.n);
            Index=[cell2mat(keys(f.c).');cell2mat(keys(g.c).')];
            Index=unique(Index,'rows');
            for i=1:(size(Index,1))
                t=Index(i,:);
                if (get(f,t)-get(g,t))~=0
                    set(h,t,get(f,t)-get(g,t));
                end
            end
        end
        %% 乘法
        function h = mtimes(f,g)
            %C = MTIMES(A,B) is called for the syntax 'A * B' when A or B is a
            if isnumeric(f)&&length(f)==1% is  scalar
                h=CZ(g.n);
                Index=cell2mat(keys(g.c).');
                for i=1:(size(Index,1))
                    t=Index(i,:);
                    set(h,t,f*get(g,t));
                end
                return
            end
            if isnumeric(g)&&length(g)==1% is  scalar
                h=g*f;
                return
            end
            N=f.n;
            h=CZ(N);
            Index1=cell2mat(keys(f.c).');
            Index2=cell2mat(keys(g.c).');
            for i=1:size(Index1,1)
                x=Index1(i,:);
%                 disp(i/size(Index1,1))
                for j=1:size(Index2,1)
                    y=Index2(j,:);
                    t=mod(x+y,N);
                    set(h,t,get(h,t)+get(f,x)*get(g,y));
                end
            end
        end
        %% 复共轭 	ctranspose(a)
        function h = ctranspose(f)
            N=f.n;
            h=CZ(N);
            Index1=cell2mat(keys(f.c).');
            for i=1:size(Index1,1)
                t=Index1(i,:);
                tp=char(N-t);
                set(h,tp,(get(f,t)'));
            end
        end
        function h=conj(f)
            h=f';
        end
        %% 乘法2 .*
        function h=times(f,g)
            %C = TIMES(A,B) is called for the syntax 'A .* B'
            %             N=f.n;
            %             m=length(N);
            %             x=sym('x',[1,m]);
            %             h= restraction(Finite_S_n(sym(f)*sym(g),x),N);
            % %             h=f*g;
            h=ProdTest(f,g);
        end
        %% Display
        function disp(f)
            S=sym(f);
            disp(S)
        end
        %% sym多项式化
        function S=sym(f)
            x=sym('x',[length(f.n),1]);
            S=sym(0);
            K=keys(f.c);
            for i=K
                ind=i{1};
                ind=double(ind);
                ind=mod(ind,f.n);
                y=get(f,ind);
                term=y*prod(x(:).^(ind(:)));
                S=S+term;
            end
        end
        
        %% Find
        function [subs,vals] = find(f)
            subs=keys(f.c);
            subs=cell2mat(subs.');
            subs=double(subs);
            vals=values(f.c);
            vals=cell2mat(vals.');
        end
        %% IFFT
        function V=CZifft(f)
            [A,B]=find(f);
            if length(f.n)>1
                F=zeros(f.n);
            else
                F=zeros(f.n,1);
                A=A+1;
                F(A+0)=B;
                V=ifft(F*length(F(:)));
                return 
            end
            N=f.n;
            for i=1:size(A,1)
                t=A(i,:);
                t=mod(t,N);
                t=t+1;
                t(t==0)=N(t==0);
                ndx=mysub2ind(f.n,t);
                F(ndx)=B(i);
            end
            V=ifftn(F*length(F(:)));
        end
        %% 在点处函数值
        function S=PointValue(this,x)
            x=x(:)';
            x=exp(2*pi*1i.*(x(:)')./(this.n));
            [a,b]=find(this);
            a=double(a);
            S=0;
            for i=1:length(b)
                t=prod(x.^[a(i,:)]);
                S=S+b(i)*t;
            end
        end
        %% Simplify
        function f=simplify(f,tol)
            %remove abs(coeffs)<tol
            if nargin==1
                tol=0;
            end
            I=keys(f.c);
            for i=1:length(I)
                if(abs(get(f,I{i}))<tol)
                    remove(f.c,I{i});
                end
            end
        end
        % min_value
        function z=min_value(f)
            V=CZifft(f);
            z=min(V(:));
        end
        
        %% isreal
        function z=isreal(f)
            [~,b]=find(f);
            z=isreal(b);
        end
    end
end




function ndx = mysub2ind(siz,v)
v1=v(1);
if length(v)>1
    v2=v(2);
    if length(v)>2
        varargin=cell(length(v)-2,1);
    end
end
for i=3:length(v)
varargin{i-2}=v(i);
end
siz = double(siz);
lensiz = length(siz);
if lensiz < 2
    error(message('MATLAB:sub2ind:InvalidSize'));
end

numOfIndInput = length(v);
if lensiz < numOfIndInput
    %Adjust for trailing singleton dimensions
    siz = [siz, ones(1,numOfIndInput-lensiz)];
elseif lensiz > numOfIndInput
    %Adjust for linear indexing on last element
    siz = [siz(1:numOfIndInput-1), prod(siz(numOfIndInput:end))];
end

if any(min(v1(:)) < 1) || any(max(v1(:)) > siz(1))
    %Verify subscripts are within range
    error(message('MATLAB:sub2ind:IndexOutOfRange'));
end

ndx = double(v1);
s = size(v1);
if numOfIndInput >= 2
    if ~isequal(s,size(v2))
        %Verify sizes of subscripts
        error(message('MATLAB:sub2ind:SubscriptVectorSize'));
    end
    if any(min(v2(:)) < 1) || any(max(v2(:)) > siz(2))
        %Verify subscripts are within range
        error(message('MATLAB:sub2ind:IndexOutOfRange'));
    end
    %Compute linear indices
    ndx = ndx + (double(v2) - 1).*siz(1);
end 
    
if numOfIndInput > 2
    %Compute linear indices
    k = cumprod(siz);
    for i = 3:numOfIndInput
        v = varargin{i-2};
        %%Input checking
        if ~isequal(s,size(v))
            %Verify sizes of subscripts
            error(message('MATLAB:sub2ind:SubscriptVectorSize'));
        end
        if (any(min(v(:)) < 1)) || (any(max(v(:)) > siz(i)))
            %Verify subscripts are within range
            error(message('MATLAB:sub2ind:IndexOutOfRange'));
        end
        ndx = ndx + (double(v)-1)*k(i-1);
    end
end
end



function v = myind2sub(siz,ndx)


nout = max(length(siz),1);
siz = double(siz);
lensiz = length(siz);

if lensiz < nout
    siz = [siz ones(1,nout-lensiz)];
elseif lensiz > nout
    siz = [siz(1:nout-1) prod(siz(nout:end))];
end

if nout > 2
    k = cumprod(siz);
    for i = nout:-1:3
        vi = rem(ndx-1, k(i-1)) + 1;
        vj = (ndx - vi)/k(i-1) + 1;
        varargout{i-2} = double(vj);
        ndx = vi;
    end
end

if nout >= 2
    vi = rem(ndx-1, siz(1)) + 1;
    v2 = double((ndx - vi)/siz(1) + 1);
    v1 = double(vi);
else
    v1 = double(ndx);
end
v=[v1,v2,cell2mat(varargout)];
end



function h=ProdTest(f,g)
N=f.n;
h=CZ(N);
[Index1,b1]=find(f);
[Index2,b2]=find(g);
y1=ones(size(Index1,1),1);
y2=ones(size(Index2,1),1);
I1=kron(Index1,y2);
I2=kron(y1,Index2);
I=mod(I1+I2,N);
I=unique(I,'rows');
I=char(I);
str=mat2cell(I,ones(1,size(I,1)),length(N));
ch=containers.Map(str,zeros(length(str),1),'uniformValues', false );
for i=1:size(Index1,1)
    x=Index1(i,:);
%     disp(i/size(Index1,1))
    t=mod(x+Index2,N);
    t=char(t);
%     str=mat2cell(t,ones(1,size(t,1)),length(N));
    coeff=b1(i)*b2;
    
    for j=1:length(coeff)
        ch(t(j,:))=ch(t(j,:))+coeff(j);
    end
end
h=setcoeff(h,ch);
end
