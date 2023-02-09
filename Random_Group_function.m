function [f,m]=Random_Group_function(N,num_terms,deg_bound)
% f is the random int-valued GP funcion
% max_{x \in G}(f(x)) <=m
if nargin==2
    deg_bound=-1;
end

if ~all(N==2)
    [f,m]=rand_int_Zn_function(N,num_terms,deg_bound);
    
else
    [f,m]=rand_int_Z2_function(N,num_terms,deg_bound);
end
V=CZifft(f);
f=f-min(V(:));
m=m-min(V(:));

end


function [f,m]=rand_int_Z2_function(N,num_terms,deg_bound)
n=length(N);
f=CZ(N);
m=0;
for i=1:num_terms
    if deg_bound~=-1
        t=zeros(1,n);
        for i=1:deg_bound
            i=randi(n);
            t(i)=t(i)+1;
        end
        for i=find(t)
            t(i)=randi(t(i));
        end
        k=randi(10);
        k=k-5;
        m=m+abs(k);
        f(t)=k;
    else
        t=rand(1,n).*N;
        t=round(t);
        t=mod(t,N);
        k=randi(10);
        k=k-5;
        m=m+abs(k);
        f(t)=k;
    end
end

end


function [f,m]=rand_int_Zn_function(N,num_terms,deg_bound)
n=length(N);
f=CZ(N);
m=0;
while(length(f))<num_terms
    t=randperm(n,randi(deg_bound));
    t=t(:)';
    Nt=N(t);
    if length(t)==1
        V=randi(10,[Nt,1]);
    else
        V=randi(10,[Nt]);
    end
    g=CZfft(V);
    m=m+max(abs(V(:)));
    [a,b]=find(g);
    for j=1:length(b)
        ai=zeros(1,n);
        ai(t)=a(j,:);
        f(ai)=f(ai)+b(j);
    end
end
end




function f=CZfft(V,tol)
if nargin==1
    tol=0;
end
if(isvector(V))
    N=length(V);
    f=CZ(N);
    F=fft(V);
    x=find(abs(F)>tol);
    for i=1:length(x)
        f(x(i)-1)=F(x(i))/N;
    end
    return
end
N=size(V);
f=CZ(N);
F=fftn(V);

for i=1:prod(N)
    if abs(F(i))>tol
        v = myind2sub(N,i);
        v=v-1;
        f(char(v))=F(i)/prod(N);
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
varargout={};
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
