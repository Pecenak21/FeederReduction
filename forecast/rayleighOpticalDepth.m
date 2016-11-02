%% Solar Resource Assessment
%  Solar Library
%
%  Title: Rayleigh Optical Depth
%  
%  Author: Bryan Urquhart
%
%
%  Description:
%    Computes the optical depth of a clean dry atmosphere
%
%
function delcda = rayleighOpticalDepth( zenith , varargin )
%% Process Input Arguments

type__ = 0;
LINKE     = type__; type__ = type__ + 1;
KASTEN    = type__; type__ = type__ + 1;
LOUCHE    = type__; type__ = type__ + 1;
MOLINEAUX = type__; %type__ = type__ + 1;

type = KASTEN;

if( ~isempty(varargin) )
  % Pass args to arg Handler
  args = argHandler(varargin);
  % Switch on args
  for i = 1:size(args,1)
    switch(lower(args{i,1}))
      case 'type'
        switch( args{i,2} )
          case 'linke'
            type = LINKE;
          case 'kasten'
            type = KASTEN;
          case 'louche'
            type = LOUCHE;
          case 'molineaux'
            type = MOLINEAUX;
          otherwise
            error('Unknown model selected.');
        end
    end
  end
end

%% Compute clear sky optical depth

am = airmass(zenith);

switch( type )
  case LINKE
    delcda = linke( am );
  case KASTEN
    delcda = kasten( am );
  case LOUCHE
    delcda = louche( am );
  case MOLINEAUX
    delcda = molineaux( am );
  otherwise
    error( 'ERROR: Unknown parameter selected!' );
end


end

%% Support functions
% Linke optical depth
function delcda = linke( am )
  delcda = 0.128 - 0.054 * log (am);
end
% Kasten optical depth
function delcda = kasten( am )
  a = [9.4 0.9];
  delcda = polyinv( a , am );
  %delcda = (9.4 + 0.9 * am).^-1;
end
% Louche optical
function delcda = louche( am )
  a = [ 6.6296, 1.7513, -0.1202, 0.0065, -0.00013 ];
  delcda = polyinv( a , am );
end
% Molineaux optical depth
function delcda = molineaux( am )
  delcda = 0.124 - 0.0656 * log(am);
end
% Inverse polynomial convenience function
function val = polyinv( a , x )
  val = 0;
  n = length(a);
  for idx = 1:n;
    val = val + a(idx) * x .^ ( idx - 1 );
  end
  val = 1./val;
end
    