
CREATE TABLE stars(
   ID INT PRIMARY KEY     NOT NULL,
   ra             REAL,
   dec		  REAL,
   catalogmag	  REAL
);

CREATE TABLE dates(
    ID INT PRIMARY KEY  NOT NULL,
    mjd      REAL,
    sunAlt   REAL,
    moonAlt  REAL,
    moonPhase REAL
);

CREATE TABLE obs(
  ID  INT PRIMARY KEY  NOT NULL,
  starID     INT,
  dateID     INT,
  alt        REAL,
  starMag    REAL,     
  starMag_err REAL,
  sky         REAL,
  filter     TEXT
);

CREATE TABLE photdiode(
 mjd    REAL,
 R      REAL,
 Y      REAL,
 Z      REAL
);

.separator ","
.import obsTable.dat obs

.import mjdTable.dat dates
.import starTable.dat stars


CREATE INDEX obsStarID on obs (starID);
CREATE INDEX obsDateID on obs (dateID);

.separator " "
.import ../photodiode/ut010215.dat photdiode
.import ../photodiode/ut012414.dat photdiode
.import ../photodiode/ut021114.dat photdiode
.import ../photodiode/ut041614.dat photdiode
.import ../photodiode/ut060614.dat photdiode
.import ../photodiode/ut091514.dat photdiode
.import ../photodiode/ut102614.dat photdiode
.import ../photodiode/ut120314.dat photdiode
.import ../photodiode/ut010315.dat photdiode
.import ../photodiode/ut012514.dat photdiode
.import ../photodiode/ut021214.dat photdiode
.import ../photodiode/ut041714.dat photdiode
.import ../photodiode/ut060714.dat photdiode
.import ../photodiode/ut091614.dat photdiode
.import ../photodiode/ut102714.dat photdiode
.import ../photodiode/ut120414.dat photdiode
.import ../photodiode/ut010415.dat photdiode
.import ../photodiode/ut012515.dat photdiode
.import ../photodiode/ut021314.dat photdiode
.import ../photodiode/ut041814.dat photdiode
.import ../photodiode/ut060814.dat photdiode
.import ../photodiode/ut091714.dat photdiode
.import ../photodiode/ut102814.dat photdiode
.import ../photodiode/ut120614.dat photdiode
.import ../photodiode/ut010515.dat photdiode
.import ../photodiode/ut012614.dat photdiode
.import ../photodiode/ut021414.dat photdiode
.import ../photodiode/ut042014.dat photdiode
.import ../photodiode/ut061214.dat photdiode
.import ../photodiode/ut092214.dat photdiode
.import ../photodiode/ut102914.dat photdiode
.import ../photodiode/ut120714.dat photdiode
.import ../photodiode/ut010615.dat photdiode
.import ../photodiode/ut012615.dat photdiode
.import ../photodiode/ut021514.dat photdiode
.import ../photodiode/ut042114.dat photdiode
.import ../photodiode/ut062614.dat photdiode
.import ../photodiode/ut092314.dat photdiode
.import ../photodiode/ut103014.dat photdiode
.import ../photodiode/ut120914.dat photdiode
.import ../photodiode/ut010915.dat photdiode
.import ../photodiode/ut012714.dat photdiode
.import ../photodiode/ut021614.dat photdiode
.import ../photodiode/ut042214.dat photdiode
.import ../photodiode/ut070314.dat photdiode
.import ../photodiode/ut092414.dat photdiode
.import ../photodiode/ut103114.dat photdiode
.import ../photodiode/ut121014.dat photdiode
.import ../photodiode/ut011015.dat photdiode
.import ../photodiode/ut012814.dat photdiode
.import ../photodiode/ut021714.dat photdiode
.import ../photodiode/ut042314.dat photdiode
.import ../photodiode/ut072114.dat photdiode
.import ../photodiode/ut092514.dat photdiode
.import ../photodiode/ut110114.dat photdiode
.import ../photodiode/ut121114.dat photdiode
.import ../photodiode/ut011115.dat photdiode
.import ../photodiode/ut012914.dat photdiode
.import ../photodiode/ut021814.dat photdiode
.import ../photodiode/ut042514.dat photdiode
.import ../photodiode/ut072214.dat photdiode
.import ../photodiode/ut092614.dat photdiode
.import ../photodiode/ut110214.dat photdiode
.import ../photodiode/ut121214.dat photdiode
.import ../photodiode/ut011214.dat photdiode
.import ../photodiode/ut012915.dat photdiode
.import ../photodiode/ut021914.dat photdiode
.import ../photodiode/ut042614.dat photdiode
.import ../photodiode/ut072614.dat photdiode
.import ../photodiode/ut092714.dat photdiode
.import ../photodiode/ut110414.dat photdiode
.import ../photodiode/ut121314.dat photdiode
.import ../photodiode/ut011215.dat photdiode
.import ../photodiode/ut013014.dat photdiode
.import ../photodiode/ut022014.dat photdiode
.import ../photodiode/ut042714.dat photdiode
.import ../photodiode/ut080114.dat photdiode
.import ../photodiode/ut092914.dat photdiode
.import ../photodiode/ut110514.dat photdiode
.import ../photodiode/ut121414.dat photdiode
.import ../photodiode/ut011314.dat photdiode
.import ../photodiode/ut013114.dat photdiode
.import ../photodiode/ut022114.dat photdiode
.import ../photodiode/ut042814.dat photdiode
.import ../photodiode/ut080214.dat photdiode
.import ../photodiode/ut093014.dat photdiode
.import ../photodiode/ut110614.dat photdiode
.import ../photodiode/ut121514.dat photdiode
.import ../photodiode/ut011315.dat photdiode
.import ../photodiode/ut013115.dat photdiode
.import ../photodiode/ut022214.dat photdiode
.import ../photodiode/ut042914.dat photdiode
.import ../photodiode/ut080314.dat photdiode
.import ../photodiode/ut100114.dat photdiode
.import ../photodiode/ut110714.dat photdiode
.import ../photodiode/ut121614.dat photdiode
.import ../photodiode/ut011414.dat photdiode
.import ../photodiode/ut020114.dat photdiode
.import ../photodiode/ut022314.dat photdiode
.import ../photodiode/ut043014.dat photdiode
.import ../photodiode/ut081814.dat photdiode
.import ../photodiode/ut100714.dat photdiode
.import ../photodiode/ut110814.dat photdiode
.import ../photodiode/ut121714.dat photdiode
.import ../photodiode/ut011415.dat photdiode
.import ../photodiode/ut020115.dat photdiode
.import ../photodiode/ut022414.dat photdiode
.import ../photodiode/ut050114.dat photdiode
.import ../photodiode/ut081914.dat photdiode
.import ../photodiode/ut100914.dat photdiode
.import ../photodiode/ut110914.dat photdiode
.import ../photodiode/ut121814.dat photdiode
.import ../photodiode/ut011514.dat photdiode
.import ../photodiode/ut020214.dat photdiode
.import ../photodiode/ut022514.dat photdiode
.import ../photodiode/ut050214.dat photdiode
.import ../photodiode/ut090114.dat photdiode
.import ../photodiode/ut101014.dat photdiode
.import ../photodiode/ut111014.dat photdiode
.import ../photodiode/ut121914.dat photdiode
.import ../photodiode/ut011515.dat photdiode
.import ../photodiode/ut020215.dat photdiode
.import ../photodiode/ut022614.dat photdiode
.import ../photodiode/ut050614.dat photdiode
.import ../photodiode/ut090214.dat photdiode
.import ../photodiode/ut101114.dat photdiode
.import ../photodiode/ut111114.dat photdiode
.import ../photodiode/ut122014.dat photdiode
.import ../photodiode/ut011614.dat photdiode
.import ../photodiode/ut020314.dat photdiode
.import ../photodiode/ut022714.dat photdiode
.import ../photodiode/ut050814.dat photdiode
.import ../photodiode/ut090314.dat photdiode
.import ../photodiode/ut101414.dat photdiode
.import ../photodiode/ut111214.dat photdiode
.import ../photodiode/ut122114.dat photdiode
.import ../photodiode/ut011615.dat photdiode
.import ../photodiode/ut020315.dat photdiode
.import ../photodiode/ut022814.dat photdiode
.import ../photodiode/ut052014.dat photdiode
.import ../photodiode/ut090414.dat photdiode
.import ../photodiode/ut101514.dat photdiode
.import ../photodiode/ut111314.dat photdiode
.import ../photodiode/ut122214.dat photdiode
.import ../photodiode/ut011714.dat photdiode
.import ../photodiode/ut020414.dat photdiode
.import ../photodiode/ut030114.dat photdiode
.import ../photodiode/ut052414.dat photdiode
.import ../photodiode/ut090514.dat photdiode
.import ../photodiode/ut101614.dat photdiode
.import ../photodiode/ut111414.dat photdiode
.import ../photodiode/ut122614.dat photdiode
.import ../photodiode/ut011715.dat photdiode
.import ../photodiode/ut020415.dat photdiode
.import ../photodiode/ut030214.dat photdiode
.import ../photodiode/ut052714.dat photdiode
.import ../photodiode/ut090614.dat photdiode
.import ../photodiode/ut101714.dat photdiode
.import ../photodiode/ut111614.dat photdiode
.import ../photodiode/ut122714.dat photdiode
.import ../photodiode/ut011814.dat photdiode
.import ../photodiode/ut020514.dat photdiode
.import ../photodiode/ut030314.dat photdiode
.import ../photodiode/ut052814.dat photdiode
.import ../photodiode/ut090714.dat photdiode
.import ../photodiode/ut101814.dat photdiode
.import ../photodiode/ut111814.dat photdiode
.import ../photodiode/ut122814.dat photdiode
.import ../photodiode/ut011815.dat photdiode
.import ../photodiode/ut020515.dat photdiode
.import ../photodiode/ut040214.dat photdiode
.import ../photodiode/ut052914.dat photdiode
.import ../photodiode/ut090814.dat photdiode
.import ../photodiode/ut101914.dat photdiode
.import ../photodiode/ut111914.dat photdiode
.import ../photodiode/ut122914.dat photdiode
.import ../photodiode/ut012014.dat photdiode
.import ../photodiode/ut020614.dat photdiode
.import ../photodiode/ut041014.dat photdiode
.import ../photodiode/ut053014.dat photdiode
.import ../photodiode/ut090914.dat photdiode
.import ../photodiode/ut102014.dat photdiode
.import ../photodiode/ut112014.dat photdiode
.import ../photodiode/ut123014.dat photdiode
.import ../photodiode/ut012015.dat photdiode
.import ../photodiode/ut020714.dat photdiode
.import ../photodiode/ut041114.dat photdiode
.import ../photodiode/ut060114.dat photdiode
.import ../photodiode/ut091014.dat photdiode
.import ../photodiode/ut102114.dat photdiode
.import ../photodiode/ut112314.dat photdiode
.import ../photodiode/ut123114.dat photdiode
.import ../photodiode/ut012114.dat photdiode
.import ../photodiode/ut020715.dat photdiode
.import ../photodiode/ut041214.dat photdiode
.import ../photodiode/ut060214.dat photdiode
.import ../photodiode/ut091114.dat photdiode
.import ../photodiode/ut102214.dat photdiode
.import ../photodiode/ut112414.dat photdiode
.import ../photodiode/ut012214.dat photdiode
.import ../photodiode/ut020814.dat photdiode
.import ../photodiode/ut041314.dat photdiode
.import ../photodiode/ut060314.dat photdiode
.import ../photodiode/ut091214.dat photdiode
.import ../photodiode/ut102314.dat photdiode
.import ../photodiode/ut112514.dat photdiode
.import ../photodiode/ut012215.dat photdiode
.import ../photodiode/ut020914.dat photdiode
.import ../photodiode/ut041414.dat photdiode
.import ../photodiode/ut060414.dat photdiode
.import ../photodiode/ut091314.dat photdiode
.import ../photodiode/ut102414.dat photdiode
.import ../photodiode/ut113014.dat photdiode

CREATE INDEX diodeDateID on photdiode (mjd);

