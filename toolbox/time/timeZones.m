function tobj = timeZones()

% timeZones    convert time to and from UTC (same as GMT or Zulu time)
%              convert time from one place to another
%   TOBJ = timeZones() creates a time zone object
%
% METHODS
%   t = timeZones()             construct time converter object
%   tz = t.zone(place)          find a time zone (-11 to +12)
%   tm = t.st2utc(time, place)  convert local standard time to UTC
%   tm = t.dst2utc(time, place) convert local daylight time to UTC
%   tm = t.utc2st(time, place)  convert utc to local standard time
%   tm = t.utc2dst(time, place) convert utc to local daylight time
%   t.places()                  all places for which time is tabulated
%
%   Most computers provide various time-related services.  Coordinated
%   Universal Time (UTC) is directly available via the Network Time
%   Protocol (NTP) (see http://www.ntp.org/). UTC is typically converted
%   to local time with a time zone offset which has to be set into the
%   computer (on WIN see the control panel Date and Time service).
%
%   When the computer is moved to a new location the setting must be
%   changed to reflect the local time. If one wants the computer to think
%   it is in Greenwich, England (that is, to show UTC), one needs to set
%   time zone 0.  If one, on the other hand, does not want to change
%   the computer setting, this application computes UTC based on a physical
%   location.  The locations are expressed as city names.  Sometimes just a
%   country, state or province name will suffice.
%   NOTE: many countries and some states sit across time zone boundaries.
%   This application will arbitrarily pick one of such ambiguous zones.
%
%   Time zone map: http://www.metoffice.gov.uk/education/images/time1pr.gif
%
% EXAMPLES
%   t = timeZones();               % construct time zone object
%   z = t.zone('Boston');          % time zone for Boston
%   assert(z == -5);
%   z = t.zone('MA');              % time zone for MA
%   assert(z == -5);
%   z = t.zone('Massachusetts');   % time zone for Massachusetts
%   assert(z == -5);
%   z = t.zone('London, England'); % time zone for London, England
%   assert(z == 0);
%   z = t.zone('London, ON');      % time zone for London, Ontario, Canada
%   assert(z == -5);
%
%   tm = now;
%   str = datestr(tm);             % text form of date/time
%   tm = t.st2utc(tm, 'Boston');   % convert Boston ST to UTC
%   tm = t.utc2st(tm, 'Boston');   % convert UTC to Boston ST
%   assert(strcmp(str, datestr(tm)));
%   str = '01-Jun-2008';           % summer day
%   tm = datenum(str);             % same units as 'now'
%   tm = t.dst2utc(tm, 'Boston');  % convert Boston DST to UTC
%   tm = t.utc2dst(tm, 'Boston');  % convert UTC to Boston DST
%   assert(strcmp(str, datestr(tm)));

%   tm = utc2st(st2utc(now, 'Boston'), 'Macau');
%   timeNowInMacau = datestr(tm);  % convert Boston time to Macau time

% http://www.infoplease.com/ipa/A0001769.html
% These values are standard local time offsets from UTC (or GMT)
persistent citiesWorld;
if isempty(citiesWorld)
	citiesWorld = {
		struct('city', 'Aberdeen, Scotland',             'time',   0)
		struct('city', 'Adelaide, Australia',            'time',   9.5)
		struct('city', 'Algiers, Algeria',               'time',   1)
		struct('city', 'Amsterdam, Netherlands',         'time',   1)
		struct('city', 'Ankara, Turkey',                 'time',   2)
		struct('city', 'Asuncion, Paraguay',             'time',  -4)
		struct('city', 'Asunci�n, Paraguay',             'time',  -4)
		struct('city', 'Athens, Greece',                 'time',   2)
		struct('city', 'Auckland, New Zealand',          'time',  12)
		struct('city', 'Bangkok, Thailand',              'time',   7)
		struct('city', 'Barcelona, Spain',               'time',   1)
		struct('city', 'Beijing, China',                 'time',   8)
		struct('city', 'Belem, Brazil',                  'time',  -3)
		struct('city', 'Bel�m, Brazil',                  'time',  -3)
		struct('city', 'Belfast, Northern Ireland',      'time',   0)
		struct('city', 'Belgrade, Serbia',               'time',   1)
		struct('city', 'Berlin, Germany',                'time',   1)
		struct('city', 'Birmingham, England',            'time',   0)
		struct('city', 'Bogota, Colombia',               'time',  -5)
		struct('city', 'Bogot�, Colombia',               'time',  -5)
		struct('city', 'Bombay, India',                  'time',   5.5)
		struct('city', 'Bordeaux, France',               'time',   1)
		struct('city', 'Bremen, Germany',                'time',   1)
		struct('city', 'Brisbane, Australia',            'time',  10)
		struct('city', 'Bristol, England',               'time',   0)
		struct('city', 'Brussels, Belgium',              'time',   1)
		struct('city', 'Bucharest, Romania',             'time',   2)
		struct('city', 'Budapest, Hungary',              'time',  -3)
		struct('city', 'Buenos Aires, Argentina',        'time',   2)
		struct('city', 'Cairo, Egypt',                   'time',   0)
		struct('city', 'Calcutta, India',                'time',   5.5)
		struct('city', 'Calgary, Canada',                'time',  -7)
		struct('city', 'Canton, China',                  'time',   8)
		struct('city', 'Cape Town, South Africa',        'time',   2)
		struct('city', 'Caracas, Venezuela',             'time',  -4)
		struct('city', 'Cayenne, French Guiana',         'time',  -3)
		struct('city', 'Chihuahua, Mexico',              'time',  -7)
		struct('city', 'Chongqing, China',               'time',   8)
		struct('city', 'Copenhagen, Denmark',            'time',   1)
		struct('city', 'Cordoba, Argentina',             'time',  -3)
		struct('city', 'C�rdoba, Argentina',             'time',  -3)
		struct('city', 'Dakar, Senegal',                 'time',   0)
		struct('city', 'Darwin, Australia',              'time',   9.5)
		struct('city', 'Djibouti, Djibouti',             'time',   3)
		struct('city', 'Dublin, Ireland',                'time',   0)
		struct('city', 'Durban, South Africa',           'time',   2)
		struct('city', 'Edinburgh, Scotland',            'time',   0)
		struct('city', 'Frankfurt, Germany',             'time',   1)
		struct('city', 'Georgetown, Guyana',             'time',  -4)
		struct('city', 'Glasgow, Scotland',              'time',   0)
		struct('city', 'Guatemala City, Guatemala',      'time',  -6)
		struct('city', 'Guayaquil, Ecuador',             'time',  -5)
		struct('city', 'Hamburg, Germany',               'time',   1)
		struct('city', 'Hammerfest, Norway',             'time',   1)
		struct('city', 'Havana, Cuba',                   'time',  -5)
		struct('city', 'Helsinki, Finland',              'time',   2)
		struct('city', 'Hobart, Tasmania',               'time',  10)
		struct('city', 'Hong Kong, China',               'time',   8)
		struct('city', 'Iquique, Chile',                 'time',  -4)
		struct('city', 'Irkutsk, Russia',                'time',   8)
		struct('city', 'Jakarta, Indonesia',             'time',   7)
		struct('city', 'Johannesburg, South Africa',     'time',   2)
		struct('city', 'Kabul, Afganistan',              'time',   4.5)
		struct('city', 'Karachi, Pakistan',              'time',   5)
		struct('city', 'Kingston, Jamaica',              'time',  -5)
		struct('city', 'Kinshasa, Congo',                'time',   1)
		struct('city', 'Kuala Lumpur, Malaysia',         'time',   8)
		struct('city', 'La Paz, Bolivia',                'time',  -4)
		struct('city', 'Leeds, England',                 'time',   0)
		struct('city', 'Lima, Peru',                     'time',  -5)
		struct('city', 'Lisbon, Portugal',               'time',   0)
		struct('city', 'Liverpool, England',             'time',   0)
		struct('city', 'London, England',                'time',   0)
		struct('city', 'London, Canada',                 'time',  -5)
		struct('city', 'Lyons, France',                  'time',   1)
		struct('city', 'Macao, China',                   'time',   8)
		struct('city', 'Macau, China',                   'time',   8)
		struct('city', 'Madrid, Spain',                  'time',   1)
		struct('city', 'Manchester, England',            'time',   0)
		struct('city', 'Manila, Philippines',            'time',   8)
		struct('city', 'Marseilles, France',             'time',   1)
		struct('city', 'Mazatlan, Mexico',               'time',  -7)
		struct('city', 'Mazatl�n, Mexico',               'time',  -7)
		struct('city', 'Mecca, Saudi Arabia',            'time',   3)
		struct('city', 'Melbourne, Australia',           'time',  10)
		struct('city', 'Mexico City, Mexico',            'time',  -6)
		struct('city', 'Milan, Italy',                   'time',   1)
		struct('city', 'Montevideo, Uruguay',            'time',  -3)
		struct('city', 'Moose Jaw, Canada',              'time',  -6)
		struct('city', 'Montreal, Canada',               'time',  -5)
		struct('city', 'Moscow, Russia',                 'time',   3)
		struct('city', 'Munich, Germany',                'time',   1)
		struct('city', 'Nagasaki, Japan',                'time',   9)
		struct('city', 'Nagoya, Japan',                  'time',   9)
		struct('city', 'Nairobi, Kenya',                 'time',   3)
		struct('city', 'Nanjing, China',                 'time',   8)
		struct('city', 'Nanking, China',                 'time',   8)
		struct('city', 'Naples, Italy',                  'time',   1)
		struct('city', 'New Delhi, India',               'time',   5.5)
		struct('city', 'Newcastle on Tyne, England',     'time',   0)
		struct('city', 'Newcastle-on-Tyne, England',     'time',   0)
		struct('city', 'Odessa, Ukraine',                'time',   2)
		struct('city', 'Osaka, Japan',                   'time',   9)
		struct('city', 'Oslo, Norway',                   'time',   1)
		struct('city', 'Ottawa, Canada',                 'time',  -5)
		struct('city', 'Panama City, Panama',            'time',  -5)
		struct('city', 'Paramaribo, Suriname',           'time',  -3)
		struct('city', 'Paris, France',                  'time',   1)
		struct('city', 'Perth, Australia',               'time',   8)
		struct('city', 'Plymouth, England',              'time',   0)
		struct('city', 'Port Moresby, Papua New Guinea', 'time',  10)
		struct('city', 'Prague, Czech Republic',         'time',   1)
		struct('city', 'Quebec, Canada',                 'time',  -5)
		struct('city', 'Rangoon, Myanmar',               'time',   6.5)
		struct('city', 'Reykjavik, Iceland',             'time',   0)
		struct('city', 'Reykjav�k, Iceland',             'time',   0)
		struct('city', 'Rio de Janeiro, Brazil',         'time',  -3)
		struct('city', 'Rome, Italy',                    'time',   1)
		struct('city', 'Salvador, Brazil',               'time',  -3)
		struct('city', 'Santiago, Chile',                'time',   8)
		struct('city', 'St Petersburg, Russia',          'time',   3)
		struct('city', 'St. Petersburg, Russia',         'time',   3)
		struct('city', 'Sao Paulo, Brazil',              'time',  -3)
		struct('city', 'S�o Paulo, Brazil',              'time',  -3)
		struct('city', 'Shanghai, China',                'time',   8)
		struct('city', 'Singapore, Singapore',           'time',   8)
		struct('city', 'Sofia, Bulgaria',                'time',   2)
		struct('city', 'Stockholm, Sweden',              'time',   1)
		struct('city', 'Sydney, Australia',              'time',  10)
		struct('city', 'Tananarive, Madagascar',         'time',   3)
		struct('city', 'Teheran, Iran',                  'time',   3.5)
		struct('city', 'Tokyo, Japan',                   'time',   9)
		struct('city', 'Toronto, Canada',                'time',  -5)
		struct('city', 'Vancouver, Canada',              'time',  -8)
		struct('city', 'Victoria, Canada',               'time',  -8)
		struct('city', 'Tripoli, Libya',                 'time',   2)
		struct('city', 'Venice, Italy',                  'time',   1)
		struct('city', 'Veracruz, Mexico',               'time',  -6)
		struct('city', 'Vienna, Austria',                'time',   1)
		struct('city', 'Vladivostok, Russia',            'time',  10)
		struct('city', 'Warsaw, Poland',                 'time',   1)
		struct('city', 'Wellington, New Zealand',        'time',  12)
		struct('city', 'Zurich, Switzerland',            'time',   1)
		struct('city', 'Z�rich, Switzerland',            'time',   1)
		};
end

% Canadian province names
persistent provinceName;
if isempty(provinceName)
	provinceName = {
		'Alberta'
		'British Columbia'
		'Manitoba'
		'New Brunswick'
		'Newfoundland and Labrador'
		'Northwest Territories'
		'Nova Scotia'
		'Nunavut'
		'Ontario'
		'Prince Edward Island'
		'Quebec'
		'Saskatchewan'
		'Yukon'
		};
end

% Canadian province abbreviations
persistent provinceAbbr;
if isempty(provinceAbbr)
	provinceAbbr = {
		'AB', 'BC', 'MB', 'NB', 'NL', 'NT', 'NS', ...
		'NU', 'ON', 'PE', 'QC', 'SK', 'YT'
		};
end

assert(numel(provinceAbbr) == numel(provinceName));  % just checking

% United States
persistent stateName;
if isempty(stateName)
	stateName ={
		'Alabama'
		'Alaska'
		'Arizona'
		'Arkansas'
		'California'
		'Colorado'
		'Connecticut'
		'Delaware'
		'Florida'
		'Georgia'
		'Hawaii'
		'Idaho'
		'Illinois'
		'Indiana'
		'Iowa'
		'Kansas'
		'Kentucky'
		'Louisiana'
		'Maine'
		'Maryland'
		'Massachusetts'
		'Michigan'
		'Minnesota'
		'Mississippi'
		'Missouri'
		'Montana'
		'Nebraska'
		'Nevada'
		'New Hampshire'
		'New Jersey'
		'New Mexico'
		'New York'
		'North Carolina'
		'North Dakota'
		'Ohio'
		'Oklahoma'
		'Oregon'
		'Pennsylvania'
		'Rhode Island'
		'South Carolina'
		'South Dakota'
		'Tennessee'
		'Texas'
		'Utah'
		'Vermont'
		'Virginia'
		'Washington'
		'West Virginia'
		'Wisconsin'
		'Wyoming'
		};
end

% State abbreviations
persistent stateAbbr;
if isempty(stateAbbr)
	stateAbbr = {
		'AL'
		'AK'
		'AZ'
		'AR'
		'CA'
		'CO'
		'CT'
		'DE'
		'FL'
		'GA'
		'HI'
		'ID'
		'IL'
		'IN'
		'IA'
		'KS'
		'KY'
		'LA'
		'ME'
		'MD'
		'MA'
		'MI'
		'MN'
		'MS'
		'MO'
		'MT'
		'NE'
		'NV'
		'NH'
		'NJ'
		'NM'
		'NY'
		'NC'
		'ND'
		'OH'
		'OK'
		'OR'
		'PA'
		'RI'
		'SC'
		'SD'
		'TN'
		'TX'
		'UT'
		'VT'
		'VA'
		'WA'
		'WV'
		'WI'
		'WY'
		};
end

assert(numel(stateAbbr) == numel(stateName));  % just checking

% http://www.infoplease.com/ipa/A0001796.html
% These values are standard local time offsets from UTC (or GMT)
persistent citiesNA;
if isempty(citiesNA)
	citiesNA = {
		struct('city', 'Albany, NY',                     'time',  -5)
		struct('city', 'Albuquerque, NM',                'time',  -7)
		struct('city', 'Amarillo, TX',                   'time',  -6)
		struct('city', 'Anchorage, AK',                  'time',  -9)
		struct('city', 'Atlanta, GA',                    'time',  -5)
		struct('city', 'Austin, TX',                     'time',  -6)
		struct('city', 'Baker, OR',                      'time',  -8)
		struct('city', 'Baltimore, MD',                  'time',  -5)
		struct('city', 'Bangor, ME',                     'time',  -5)
		struct('city', 'Birmingham, AL',                 'time',  -6)
		struct('city', 'Bismarck, ND',                   'time',  -6)
		struct('city', 'Boise, ID',                      'time',  -7)
		struct('city', 'Boston, MA',                     'time',  -5)
		struct('city', 'Buffalo, NY',                    'time',  -5)
		struct('city', 'Calgary, AB',                    'time',  -7)
		struct('city', 'Carlsbad, NM',                   'time',  -7)
		struct('city', 'Charleston, SC',                 'time',  -5)
		struct('city', 'Charleston, WV',                 'time',  -5)
		struct('city', 'Charlotte, NC',                  'time',  -5)
		struct('city', 'Charlottetown, PE',              'time',  -4)
		struct('city', 'Chesterfield Inlet, NU',         'time',  -6)
		struct('city', 'Cheyenne, WY',                   'time',  -7)
		struct('city', 'Chicago, IL',                    'time',  -6)
		struct('city', 'Cincinnati, OH',                 'time',  -5)
		struct('city', 'Cleveland, OH',                  'time',  -5)
		struct('city', 'Columbia, SC',                   'time',  -5)
		struct('city', 'Columbus, OH',                   'time',  -5)
		struct('city', 'Dallas, TX',                     'time',  -6)
		struct('city', 'Denver, CO',                     'time',  -7)
		struct('city', 'Des Moines, IA',                 'time',  -6)
		struct('city', 'Detroit, MI',                    'time',  -5)
		struct('city', 'Dover, DE',                      'time',  -5)
		struct('city', 'Dubuque, IA',                    'time',  -6)
		struct('city', 'Duluth, MN',                     'time',  -6)
		struct('city', 'Eastport, ME',                   'time',  -5)
		struct('city', 'Edmonton, AB',                   'time',  -7)
		struct('city', 'El Centro, CA',                  'time',  -8)
		struct('city', 'El Paso, TX',                    'time',  -7)
		struct('city', 'Eugene, OR',                     'time',  -8)
		struct('city', 'Fargo, ND',                      'time',  -6)
		struct('city', 'Flagstaff, AZ',                  'time',  -7)
		struct('city', 'Fort Worth, TX',                 'time',  -6)
		struct('city', 'Fresno, CA',                     'time',  -8)
		struct('city', 'Grand Junction, CO',             'time',  -7)
		struct('city', 'Grand Rapids, MI',               'time',  -5)
		struct('city', 'Halifax, NS',                    'time',  -4)
		struct('city', 'Havre, MT',                      'time',  -7)
		struct('city', 'Helena, MT',                     'time',  -7)
		struct('city', 'Honolulu, HI',                   'time',  -10)
		struct('city', 'Hot Springs, AR',                'time',  -6)
		struct('city', 'Houston, TX',                    'time',  -6)
		struct('city', 'Idaho Falls, ID',                'time',  -7)
		struct('city', 'Indianapolis, IN',               'time',  -5)
		struct('city', 'Jackson, MS',                    'time',  -6)
		struct('city', 'Jacksonville, FL',               'time',  -5)
		struct('city', 'Juneau, AK',                     'time',  -9)
		struct('city', 'Kansas City, MO',                'time',  -6)
		struct('city', 'Key West, FL',                   'time',  -5)
		struct('city', 'Kingston, ON',                   'time',  -5)
		struct('city', 'Klamath Falls, OR',              'time',  -8)
		struct('city', 'Knoxville, TN',                  'time',  -5)
		struct('city', 'Las Vegas, NV',                  'time',  -8)
		struct('city', 'Lewiston, ID',                   'time',  -8)
		struct('city', 'Lincoln, NB',                    'time',  -6)
		struct('city', 'London, ON',                     'time',  -5)
		struct('city', 'Long Beach, CA',                 'time',  -8)
		struct('city', 'Los Angeles, CA',                'time',  -8)
		struct('city', 'Louisville, KY',                 'time',  -5)
		struct('city', 'Madison, WI',                    'time',  -6)
		struct('city', 'Manchester, NH',                 'time',  -5)
		struct('city', 'Memphis, TN',                    'time',  -6)
		struct('city', 'Miami, FL',                      'time',  -5)
		struct('city', 'Milwaukee, WS',                  'time',  -6)
		struct('city', 'Minneapolis, MN',                'time',  -6)
		struct('city', 'Mobile, AL',                     'time',  -6)
		struct('city', 'Montgomery, AL',                 'time',  -6)
		struct('city', 'Montpelier, VT',                 'time',  -5)
		struct('city', 'Montreal, QC',                   'time',  -5)
		struct('city', 'Moose Jaw, SK',                  'time',  -6)
		struct('city', 'Nashville, TN',                  'time',  -6)
		struct('city', 'Nelson, BC',                     'time',  -8)
		struct('city', 'Newark, NJ',                     'time',  -5)
		struct('city', 'New Haven, CT',                  'time',  -5)
		struct('city', 'New Orleans, LA',                'time',  -6)
		struct('city', 'New York, NY',                   'time',  -5)
		struct('city', 'Nome, AK',                       'time',  -9)
		struct('city', 'Oakland, CA',                    'time',  -8)
		struct('city', 'Oklahoma City, OK',              'time',  -6)
		struct('city', 'Omaha, NE',                      'time',  -6)
		struct('city', 'Ottawa, ON',                     'time',  -5)
		struct('city', 'Philadelphia, PA',               'time',  -5)
		struct('city', 'Phoenix, AZ',                    'time',  -7)
		struct('city', 'Pierre, SD',                     'time',  -6)
		struct('city', 'Pittsburgh, PA',                 'time',  -5)
		struct('city', 'Portland, ME',                   'time',  -5)
		struct('city', 'Portland, OR',                   'time',  -8)
		struct('city', 'Providence, RI',                 'time',  -5)
		struct('city', 'Quebec, QC',                     'time',  -5)
		struct('city', 'Raleigh, NC',                    'time',  -5)
		struct('city', 'Reno, NV',                       'time',  -8)
		struct('city', 'Richfield, UT',                  'time',  -7)
		struct('city', 'Richmond, VA',                   'time',  -5)
		struct('city', 'Roanoke, VA',                    'time',  -5)
		struct('city', 'Sacramento, CA',                 'time',  -8)
		struct('city', 'St. John''s, NL',                'time',  -3.5)
		struct('city', 'St Louis, MO',                   'time',  -6)
		struct('city', 'Salt Lake City, UT',             'time',  -7)
		struct('city', 'San Antonio, TX',                'time',  -6)
		struct('city', 'San Diego, CA',                  'time',  -8)
		struct('city', 'San Francisco, CA',              'time',  -8)
		struct('city', 'San Jose, CA',                   'time',  -8)
		struct('city', 'San Juan, PR',                   'time',  -4)
		struct('city', 'Santa Fe, NM',                   'time',  -7)
		struct('city', 'Savannah, GA',                   'time',  -5)
		struct('city', 'Seattle, WA',                    'time',  -8)
		struct('city', 'Shreveport, LA',                 'time',  -6)
		struct('city', 'Sioux Falls, SD',                'time',  -6)
		struct('city', 'Sitka, AK',                      'time',  -9)
		struct('city', 'Spokane, WA',                    'time',  -8)
		struct('city', 'Springfield, IL',                'time',  -6)
		struct('city', 'Springfield, MA',                'time',  -5)
		struct('city', 'Springfield, MO',                'time',  -6)
		struct('city', 'Syracuse, NY',                   'time',  -5)
		struct('city', 'Tampa, FL',                      'time',  -5)
		struct('city', 'Toledo, OH',                     'time',  -5)
		struct('city', 'Toronto, ON',                    'time',  -5)
		struct('city', 'Tulsa, OK',                      'time',  -6)
		struct('city', 'Vancouver, BC',                  'time',  -8)
		struct('city', 'Victoria, BC',                   'time',  -8)
		struct('city', 'Virginia Beach, VA',             'time',  -5)
		struct('city', 'Washington, DC',                 'time',  -5)
		struct('city', 'Whitehorse, YT',                 'time',  -8)
		struct('city', 'Wichita, KS',                    'time',  -6)
		struct('city', 'Wilmington, NC',                 'time',  -5)
		struct('city', 'Winnipeg, MB',                   'time',  -6)
		struct('city', 'Yellowknife, NT',                'time',  -8)
		};
end

% set up memoization for repeated calls
memo = ' -- not a place -- ';
memozone = 0;

tobj = public();                               % prepare public interface
%return;

% ------------------------------ methods ----------------------

	function tz = zone(place)
		if strcmp(place, memo)                       % use previous value
			tz = memozone;
			return
		end
		
		% table look up
		np = numel(place);                           % num chars
		for jj=1:numel(stateName)                    % use state abbreviation
			if strcmp(place, stateName{jj})            % in place of state name
				tz = zone(stateAbbr{jj});
				return;
			end
		end
		for jj=1:numel(provinceName)                 % use province abbr
			if strcmp(place, provinceName{jj})         % in place of name
				tz = zone(provinceAbbr{jj});
				return;
			end
		end
		
		for jj=1:numel(citiesNA)                     % find North American city
			entry = citiesNA{jj}.city;                 % table entry
			ne = numel(entry);                         % num chars
			if ne >= np                                % skip short entries
				candidates = strfind(entry, place);      % matches
				for kk = candidates                      % check each
					if np == ne                            ... exact match
							|| np == 2 && strcmp(entry(end-1:end), place) ... state code
							|| ne > np && kk == 1 && entry(np+1) == ','  % city name match
						tz = citiesNA{jj}.time;              % get zone
						memo = place; memozone = tz;         % in case of repeated call
						return;
					end
				end
			end
		end
		
		for jj=1:numel(citiesWorld)
			entry = citiesWorld{jj}.city;
			ne = numel(entry);                         % num chars
			if ne >= np                                % skip short entries
				candidates = strfind(entry, place);      % matches
				for kk = candidates                      % check each
					if np == ne                            ... exact match
							|| ne > np && kk == 1 && entry(np+1) == ','  ... city name match
							|| entry(kk-1) == ' ' && strcmp(place, entry(end-np+1:end))
						tz = citiesWorld{jj}.time;           % get zone
						memo = place; memozone = tz;         % in case of repeated call
						return;
					end
				end
			end
		end
		error('''%s'' is not in time zone tables', place);
	end

% convert local daylight saving time to UTC
	function utc = dst2utc(time, place)
		utc = time-1/24-zone(place)/24;
	end

% convert local standard time to UTC
	function utc = st2utc(time, place)
		utc = time-zone(place)/24;
	end

% convert UTC to local daylight saving time
	function utc = utc2dst(time, place)
		utc = time+1/24+zone(place)/24;
	end

% convert UTC to local standard time
	function utc = utc2st(time, place)
		utc = time+zone(place)/24;
	end

	function places()
		for jj = 1:numel(citiesNA)
			fprintf('%s\n', citiesNA{jj}.city);
		end
		fprintf('\n');
		for jj = 1:numel(citiesWorld)
			fprintf('%s\n', citiesWorld{jj}.city);
		end
		fprintf('\n');
		for jj = 1:numel(stateName)
			fprintf('%s\n', stateName{jj});
		end
		fprintf('\n');
		for jj = 1:numel(stateAbbr)
			fprintf('%s\n', stateAbbr{jj});
		end
		fprintf('\n');
		for jj = 1:numel(provinceName)
			fprintf('%s\n', provinceName{jj});
		end
		fprintf('\n');
		for jj = 1:numel(provinceName)
			fprintf('%s\n', provinceAbbr{jj});
		end
		fprintf('\n');
	end

% assemble the public interface
	function o = public()
		o = struct(...
			'zone',    @zone,...
			'st2utc',  @st2utc,...
			'dst2utc', @dst2utc,...
			'utc2st',  @utc2st,...
			'utc2dst', @utc2dst,...
			'places',  @places...
			);
	end

end