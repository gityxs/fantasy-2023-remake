(function (global, factory) {
    typeof exports === 'object' && typeof module !== 'undefined' ? module.exports = factory() :
    typeof define === 'function' && define.amd ? define(factory) :
    (global = global || self, global.Decimal = factory());
  }(this, function () { 'use strict';
  
    var padEnd = function (string, maxLength, fillString) {
  
      if (string === null || maxLength === null) {
        return string;
      }
  
      var result    = String(string);
      var targetLen = typeof maxLength === 'number'
        ? maxLength
        : parseInt(maxLength, 10);
  
      if (isNaN(targetLen) || !isFinite(targetLen)) {
        return result;
      }
  
  
      var length = result.length;
      if (length >= targetLen) {
        return result;
      }
  
  
      var filled = fillString === null ? '' : String(fillString);
      if (filled === '') {
        filled = ' ';
      }
  
  
      var fillLen = targetLen - length;
  
      while (filled.length < fillLen) {
        filled += filled;
      }
  
      var truncated = filled.length > fillLen ? filled.substr(0, fillLen) : filled;
  
      return result + truncated;
    };
  
    var MAX_SIGNIFICANT_DIGITS = 17; //Maximum number of digits of precision to assume in Number
  
    var EXP_LIMIT = 9e15; //If we're ABOVE this value, increase a layer. (9e15 is close to the largest integer that can fit in a Number.)
    
    var LAYER_DOWN = Math.log10(9e15); //If we're BELOW this value, drop down a layer. About 15.954.
    
    var FIRST_NEG_LAYER = 1/9e15; //At layer 0, smaller non-zero numbers than this become layer 1 numbers with negative mag. After that the pattern continues as normal.
  
    var NUMBER_EXP_MAX = 308; //The largest exponent that can appear in a Number, though not all mantissas are valid here.
  
    var NUMBER_EXP_MIN = -324; //The smallest exponent that can appear in a Number, though not all mantissas are valid here.
    
    var MAX_ES_IN_A_ROW = 5; //For default toString behaviour, when to swap from eee... to (e^n) syntax.

    var powerOf10 = function () {
      // We need this lookup table because Math.pow(10, exponent)
      // when exponent's absolute value is large is slightly inaccurate.
      // You can fix it with the power of math... or just make a lookup table.
      // Faster AND simpler
      var powersOf10 = [];
  
      for (var i = NUMBER_EXP_MIN + 1; i <= NUMBER_EXP_MAX; i++) {
        powersOf10.push(Number("1e" + i));
      }
  
      var indexOf0InPowersOf10 = 323;
      return function (power) {
        return powersOf10[power + indexOf0InPowersOf10];
      };
    }();
  
    var D = function D(value) {
      return Decimal.fromValue_noAlloc(value);
    };
  
    var FC = function FC(sign, layer, mag) {
      return Decimal.fromComponents(sign, layer, mag);
    };
  
    var FC_NN = function FC_NN(sign, layer, mag) {
      return Decimal.fromComponents_noNormalize(sign, layer, mag);
    };
    
    var ME = function ME(mantissa, exponent) {
      return Decimal.fromMantissaExponent(mantissa, exponent);
    };
  
    var ME_NN = function ME_NN(mantissa, exponent) {
      return Decimal.fromMantissaExponent_noNormalize(mantissa, exponent);
    };
    
    var decimalPlaces = function decimalPlaces(value, places) {
      var len = places + 1;
      var numDigits = Math.ceil(Math.log10(Math.abs(value)));
      var rounded = Math.round(value * Math.pow(10, len - numDigits)) * Math.pow(10, numDigits - len);
      return parseFloat(rounded.toFixed(Math.max(len - numDigits, 0)));
    };
    
    var f_maglog10 = function(n) {
      return Math.sign(n)*Math.log10(Math.abs(n));
    }
    
    //from HyperCalc source code
    var f_gamma = function(n) {
      if (!isFinite(n)) { return n; }
      if (n < -50)
      {
        if (n === Math.trunc(n)) { return Number.NEGATIVE_INFINITY; }
        return 0;
      }
      
      var scal1 = 1;
      while (n < 10)
      {
        scal1 = scal1*n;
        ++n;
      }
      
      n -= 1;
      var l = 0.9189385332046727; //0.5*Math.log(2*Math.PI)
      l = l + (n+0.5)*Math.log(n);
      l = l - n;
      var n2 = n*n;
      var np = n;
      l = l+1/(12*np);
      np = np*n2;
      l = l+1/(360*np);
      np = np*n2;
      l = l+1/(1260*np);
      np = np*n2;
      l = l+1/(1680*np);
      np = np*n2;
      l = l+1/(1188*np);
      np = np*n2;
      l = l+691/(360360*np);
      np = np*n2;
      l = l+7/(1092*np);
      np = np*n2;
      l = l+3617/(122400*np);
  
      return Math.exp(l)/scal1;
    };
    
    var twopi = 6.2831853071795864769252842;  // 2*pi
    var EXPN1 = 0.36787944117144232159553;  // exp(-1)
    var OMEGA = 0.56714329040978387299997;  // W(1, 0)
    //from https://math.stackexchange.com/a/465183
    // The evaluation can become inaccurate very close to the branch point
    var f_lambertw = function(z, tol = 1e-10) {
      var w;
      var wn;
  
      if (!Number.isFinite(z)) { return z; }
      if (z === 0)
      {
        return z;
      }
      if (z === 1)
      {
        return OMEGA;
      }
  
      if (z < 10)
      {
        w = 0;
      }
      else
      {
        w = Math.log(z)-Math.log(Math.log(z));
      }
  
      for (var i = 0; i < 100; ++i)
      {
        wn = (z * Math.exp(-w) + w * w)/(w + 1);
        if (Math.abs(wn - w) < tol*Math.abs(wn))
        {
          return wn;
        }
        else
        {
          w = wn;
        }
      }
  
      throw Error("Iteration failed to converge: " + z);
      //return Number.NaN;
    }
    
    var Decimal =
    /** @class */
    function () {
    
      function Decimal(value) {
        
        this.sign = Number.NaN;
        this.layer = Number.NaN;
        this.mag = Number.NaN;
  
        if (value instanceof Decimal) {
          this.fromDecimal(value);
        } else if (typeof value === "number") {
          this.fromNumber(value);
        } else if (typeof value === "string") {
          this.fromString(value);
        } else {
          this.sign = 0;
          this.layer = 0;
          this.mag = 0;
        }
      }
  
      Object.defineProperty(Decimal.prototype, "m", {
        get: function get() {
          if (this.sign === 0)
          {
            return 0;
          }
          else if (this.layer === 0)
          {
            var exp = Math.floor(Math.log10(this.mag));
            //handle special case 5e-324
            var man;
            if (this.mag === 5e-324)
            {
              man = 5;
            }
            else
            {
              man = this.mag / powerOf10(exp);
            }
            return this.sign*man;
          }
          else if (this.layer === 1)
          {
            var residue = this.mag-Math.floor(this.mag);
            return this.sign*Math.pow(10, residue);
          }
          else
          {
            //mantissa stops being relevant past 1e9e15 / ee15.954
            return this.sign;
          }
        },
        set: function set(value) {
          if (this.layer <= 2)
          {
            this.fromMantissaExponent(value, this.e);
          }
          else
          {
            //don't even pretend mantissa is meaningful
            this.sign = Math.sign(value);
            if (this.sign === 0) { this.layer === 0; this.exponent === 0; }
          }
        },
        enumerable: true,
        configurable: true
      });
      Object.defineProperty(Decimal.prototype, "e", {
        get: function get() {
          if (this.sign === 0)
          {
            return 0;
          }
          else if (this.layer === 0)
          {
            return Math.floor(Math.log10(this.mag));
          }
          else if (this.layer === 1)
          {
            return Math.floor(this.mag);
          }
          else if (this.layer === 2)
          {
            return Math.floor(Math.sign(this.mag)*Math.pow(10, Math.abs(this.mag)));
          }
          else
          {
            return this.mag*Number.POSITIVE_INFINITY;
          }
        },
        set: function set(value) {
          this.fromMantissaExponent(this.m, value);
        },
        enumerable: true,
        configurable: true
      });
      Object.defineProperty(Decimal.prototype, "s", {
        get: function get() {
          return this.sign;
        },
        set: function set(value) {
          if (value === 0) {
            this.sign = 0;
            this.layer = 0;
            this.mag = 0;
          }
          else
          {
            this.sign = value;
          }
        },
        enumerable: true,
        configurable: true
      });
      
      Object.defineProperty(Decimal.prototype, "mantissa", {
        get: function get() {
          return this.m;
        },
        set: function set(value) {
          this.m = value;
        },
        enumerable: true,
        configurable: true
      });
  
      Object.defineProperty(Decimal.prototype, "exponent", {
        get: function get() {
          return this.e;
        },
        set: function set(value) {
          this.e = value;
        },
        enumerable: true,
        configurable: true
      });
  
      Decimal.fromComponents = function (sign, layer, mag) {
        return new Decimal().fromComponents(sign, layer, mag);
      };
  
      Decimal.fromComponents_noNormalize = function (sign, layer, mag) {
        return new Decimal().fromComponents_noNormalize(sign, layer, mag);
      };
      
      Decimal.fromMantissaExponent = function (mantissa, exponent) {
        return new Decimal().fromMantissaExponent(mantissa, exponent);
      };
  
      Decimal.fromMantissaExponent_noNormalize = function (mantissa, exponent) {
        return new Decimal().fromMantissaExponent_noNormalize(mantissa, exponent);
      };
      
      Decimal.fromDecimal = function (value) {
        return new Decimal().fromDecimal(value);
      };
  
      Decimal.fromNumber = function (value) {
        return new Decimal().fromNumber(value);
      };
  
      Decimal.fromString = function (value) {
        return new Decimal().fromString(value);
      };
  
      Decimal.fromValue = function (value) {
        return new Decimal().fromValue(value);
      };
  
      Decimal.fromValue_noAlloc = function (value) {
        return value instanceof Decimal ? value : new Decimal(value);
      };
      
      Decimal.abs = function (value) {
        return D(value).abs();
      };
  
      Decimal.neg = function (value) {
        return D(value).neg();
      };
  
      Decimal.negate = function (value) {
        return D(value).neg();
      };
  
      Decimal.negated = function (value) {
        return D(value).neg();
      };
  
      Decimal.sign = function (value) {
        return D(value).sign();
      };
  
      Decimal.sgn = function (value) {
        return D(value).sign();
      };
  
      Decimal.round = function (value) {
        return D(value).round();
      };
  
      Decimal.floor = function (value) {
        return D(value).floor();
      };
  
      Decimal.ceil = function (value) {
        return D(value).ceil();
      };
  
      Decimal.trunc = function (value) {
        return D(value).trunc();
      };
  
      Decimal.add = function (value, other) {
        return D(value).add(other);
      };
  
      Decimal.plus = function (value, other) {
        return D(value).add(other);
      };
  
      Decimal.sub = function (value, other) {
        return D(value).sub(other);
      };
  
      Decimal.subtract = function (value, other) {
        return D(value).sub(other);
      };
  
      Decimal.minus = function (value, other) {
        return D(value).sub(other);
      };
  
      Decimal.mul = function (value, other) {
        return D(value).mul(other);
      };
  
      Decimal.multiply = function (value, other) {
        return D(value).mul(other);
      };
  
      Decimal.times = function (value, other) {
        return D(value).mul(other);
      };
  
      Decimal.div = function (value, other) {
        return D(value).div(other);
      };
  
      Decimal.divide = function (value, other) {
        return D(value).div(other);
      };
  
      Decimal.recip = function (value) {
        return D(value).recip();
      };
  
      Decimal.reciprocal = function (value) {
        return D(value).recip();
      };
  
      Decimal.reciprocate = function (value) {
        return D(value).reciprocate();
      };
  
      Decimal.cmp = function (value, other) {
        return D(value).cmp(other);
      };
  
      Decimal.cmpabs = function (value, other) {
        return D(value).cmpabs(other);
      };
      
      Decimal.compare = function (value, other) {
        return D(value).cmp(other);
      };
  
      Decimal.eq = function (value, other) {
        return D(value).eq(other);
      };
  
      Decimal.equals = function (value, other) {
        return D(value).eq(other);
      };
  
      Decimal.neq = function (value, other) {
        return D(value).neq(other);
      };
  
      Decimal.notEquals = function (value, other) {
        return D(value).notEquals(other);
      };
  
      Decimal.lt = function (value, other) {
        return D(value).lt(other);
      };
  
      Decimal.lte = function (value, other) {
        return D(value).lte(other);
      };
  
      Decimal.gt = function (value, other) {
        return D(value).gt(other);
      };
  
      Decimal.gte = function (value, other) {
        return D(value).gte(other);
      };
  
      Decimal.max = function (value, other) {
        return D(value).max(other);
      };
      
      Decimal.min = function (value, other) {
        return D(value).min(other);
      };
  
      Decimal.minabs = function (value, other) {
        return D(value).minabs(other);
      };
      
      Decimal.maxabs = function (value, other) {
        return D(value).maxabs(other);
      };
      
      Decimal.clamp = function(value, min, max) {
        return D(value).clamp(min, max);
      }
      
      Decimal.clampMin = function(value, min) {
        return D(value).clampMin(min);
      }
      
      Decimal.clampMax = function(value, max) {
        return D(value).clampMax(max);
      }
  
      Decimal.cmp_tolerance = function (value, other, tolerance) {
        return D(value).cmp_tolerance(other, tolerance);
      };
  
      Decimal.compare_tolerance = function (value, other, tolerance) {
        return D(value).cmp_tolerance(other, tolerance);
      };
  
      Decimal.eq_tolerance = function (value, other, tolerance) {
        return D(value).eq_tolerance(other, tolerance);
      };
  
      Decimal.equals_tolerance = function (value, other, tolerance) {
        return D(value).eq_tolerance(other, tolerance);
      };
  
      Decimal.neq_tolerance = function (value, other, tolerance) {
        return D(value).neq_tolerance(other, tolerance);
      };
  
      Decimal.notEquals_tolerance = function (value, other, tolerance) {
        return D(value).notEquals_tolerance(other, tolerance);
      };
  
      Decimal.lt_tolerance = function (value, other, tolerance) {
        return D(value).lt_tolerance(other, tolerance);
      };
  
      Decimal.lte_tolerance = function (value, other, tolerance) {
        return D(value).lte_tolerance(other, tolerance);
      };
  
      Decimal.gt_tolerance = function (value, other, tolerance) {
        return D(value).gt_tolerance(other, tolerance);
      };
  
      Decimal.gte_tolerance = function (value, other, tolerance) {
        return D(value).gte_tolerance(other, tolerance);
      };
  
      Decimal.pLog10 = function (value) {
        return D(value).pLog10();
      };
      
      Decimal.absLog10 = function (value) {
        return D(value).absLog10();
      };
      
      Decimal.log10 = function (value) {
        return D(value).log10();
      };
  
      Decimal.log = function (value, base) {
        return D(value).log(base);
      };
  
      Decimal.log2 = function (value) {
        return D(value).log2();
      };
  
      Decimal.ln = function (value) {
        return D(value).ln();
      };
  
      Decimal.logarithm = function (value, base) {
        return D(value).logarithm(base);
      };
  
      Decimal.pow = function (value, other) {
        return D(value).pow(other);
      };
      
      Decimal.pow10 = function (value) {
        return D(value).pow10();
      };
      
      Decimal.root = function (value, other) {
        return D(value).root(other);
      };
      
      Decimal.factorial = function (value, other) {
        return D(value).factorial();
      };
      
      Decimal.gamma = function (value, other) {
        return D(value).gamma();
      };
      
      Decimal.lngamma = function (value, other) {
        return D(value).lngamma();
      };
  
      Decimal.exp = function (value) {
        return D(value).exp();
      };
  
      Decimal.sqr = function (value) {
        return D(value).sqr();
      };
  
      Decimal.sqrt = function (value) {
        return D(value).sqrt();
      };
  
      Decimal.cube = function (value) {
        return D(value).cube();
      };
  
      Decimal.cbrt = function (value) {
        return D(value).cbrt();
      };
      
      Decimal.tetrate = function (value, height = 2, payload = FC_NN(1, 0, 1)) {
        return D(value).tetrate(height, payload);
      }
      
      Decimal.iteratedexp = function (value, height = 2, payload = FC_NN(1, 0, 1)) {
        return D(value).iteratedexp(height, payload);
      }
      
      Decimal.iteratedlog = function (value, base = 10, times = 1) {
        return D(value).iteratedlog(base, times);
      }
      
      Decimal.layeradd10 = function (value, diff) {
        return D(value).layeradd10(diff);
      }
      
       Decimal.layeradd = function (value, diff, base = 10) {
        return D(value).layeradd(diff, base);
      }
      
      Decimal.slog = function (value, base = 10) {
        return D(value).slog(base);
      }
      
      Decimal.lambertw = function(value) {
        return D(value).lambertw();
      }
      
      Decimal.ssqrt = function(value) {
        return D(value).ssqrt();
      }
      
      Decimal.pentate = function (value, height = 2, payload = FC_NN(1, 0, 1)) {
        return D(value).pentate(height, payload);
      }
      
      /**
       * If you're willing to spend 'resourcesAvailable' and want to buy something
       * with exponentially increasing cost each purchase (start at priceStart,
       * multiply by priceRatio, already own currentOwned), how much of it can you buy?
       * Adapted from Trimps source code.
       */
  
  
      Decimal.affordGeometricSeries = function (resourcesAvailable, priceStart, priceRatio, currentOwned) {
        return this.affordGeometricSeries_core(D(resourcesAvailable), D(priceStart), D(priceRatio), currentOwned);
      };
      /**
       * How much resource would it cost to buy (numItems) items if you already have currentOwned,
       * the initial price is priceStart and it multiplies by priceRatio each purchase?
       */
  
  
      Decimal.sumGeometricSeries = function (numItems, priceStart, priceRatio, currentOwned) {
        return this.sumGeometricSeries_core(numItems, D(priceStart), D(priceRatio), currentOwned);
      };
      /**
       * If you're willing to spend 'resourcesAvailable' and want to buy something with additively
       * increasing cost each purchase (start at priceStart, add by priceAdd, already own currentOwned),
       * how much of it can you buy?
       */
  
  
      Decimal.affordArithmeticSeries = function (resourcesAvailable, priceStart, priceAdd, currentOwned) {
        return this.affordArithmeticSeries_core(D(resourcesAvailable), D(priceStart), D(priceAdd), D(currentOwned));
      };
      /**
       * How much resource would it cost to buy (numItems) items if you already have currentOwned,
       * the initial price is priceStart and it adds priceAdd each purchase?
       * Adapted from http://www.mathwords.com/a/arithmetic_series.htm
       */
  
  
      Decimal.sumArithmeticSeries = function (numItems, priceStart, priceAdd, currentOwned) {
        return this.sumArithmeticSeries_core(D(numItems), D(priceStart), D(priceAdd), D(currentOwned));
      };
      /**
       * When comparing two purchases that cost (resource) and increase your resource/sec by (deltaRpS),
       * the lowest efficiency score is the better one to purchase.
       * From Frozen Cookies:
       * http://cookieclicker.wikia.com/wiki/Frozen_Cookies_(JavaScript_Add-on)#Efficiency.3F_What.27s_that.3F
       */
  
  
      Decimal.efficiencyOfPurchase = function (cost, currentRpS, deltaRpS) {
        return this.efficiencyOfPurchase_core(D(cost), D(currentRpS), D(deltaRpS));
      };
  
      Decimal.randomDecimalForTesting = function (maxLayers) {
        // NOTE: This doesn't follow any kind of sane random distribution, so use this for testing purposes only.
        //5% of the time, return 0
        if (Math.random() * 20 < 1) {
          return FC_NN(0, 0, 0);
        }
        
        var randomsign = Math.random() > 0.5 ? 1 : -1;
        
        //5% of the time, return 1 or -1
        if (Math.random() * 20 < 1) {
          return FC_NN(randomsign, 0, 1);
        }
        
        //pick a random layer
        var layer = Math.floor(Math.random()*(maxLayers+1));
  
        var randomexp = layer === 0 ? Math.random()*616-308 : Math.random()*16;
        //10% of the time, make it a simple power of 10
        if (Math.random() > 0.9) { randomexp = Math.trunc(randomexp); }
        var randommag = Math.pow(10, randomexp);
        //10% of the time, trunc mag
        if (Math.random() > 0.9) { randommag = Math.trunc(randommag); }
        return FC(randomsign, layer, randommag);
      };
  
      Decimal.affordGeometricSeries_core = function (resourcesAvailable, priceStart, priceRatio, currentOwned) {
        var actualStart = priceStart.mul(priceRatio.pow(currentOwned));
        return Decimal.floor(resourcesAvailable.div(actualStart).mul(priceRatio.sub(1)).add(1).log10().div(priceRatio.log10()));
      };
  
      Decimal.sumGeometricSeries_core = function (numItems, priceStart, priceRatio, currentOwned) {
        return priceStart.mul(priceRatio.pow(currentOwned)).mul(Decimal.sub(1, priceRatio.pow(numItems))).div(Decimal.sub(1, priceRatio));
      };
  
      Decimal.affordArithmeticSeries_core = function (resourcesAvailable, priceStart, priceAdd, currentOwned) {
        // n = (-(a-d/2) + sqrt((a-d/2)^2+2dS))/d
        // where a is actualStart, d is priceAdd and S is resourcesAvailable
        // then floor it and you're done!
        var actualStart = priceStart.add(currentOwned.mul(priceAdd));
        var b = actualStart.sub(priceAdd.div(2));
        var b2 = b.pow(2);
        return b.neg().add(b2.add(priceAdd.mul(resourcesAvailable).mul(2)).sqrt()).div(priceAdd).floor();
      };
  
      Decimal.sumArithmeticSeries_core = function (numItems, priceStart, priceAdd, currentOwned) {
        var actualStart = priceStart.add(currentOwned.mul(priceAdd)); // (n/2)*(2*a+(n-1)*d)
  
        return numItems.div(2).mul(actualStart.mul(2).plus(numItems.sub(1).mul(priceAdd)));
      };
  
      Decimal.efficiencyOfPurchase_core = function (cost, currentRpS, deltaRpS) {
        return cost.div(currentRpS).add(cost.div(deltaRpS));
      };
      
      Decimal.prototype.normalize = function () {
        /*
        PSEUDOCODE:
        Whenever we are partially 0 (sign is 0 or mag and layer is 0), make it fully 0.
        Whenever we are at or hit layer 0, extract sign from negative mag.
        If layer === 0 and mag < FIRST_NEG_LAYER (1/9e15), shift to 'first negative layer' (add layer, log10 mag).
        While abs(mag) > EXP_LIMIT (9e15), layer += 1, mag = maglog10(mag).
        While abs(mag) < LAYER_DOWN (15.954) and layer > 0, layer -= 1, mag = pow(10, mag).
        
        When we're done, all of the following should be true OR one of the numbers is not IsFinite OR layer is not IsInteger (error state):
        Any 0 is totally zero (0, 0, 0).
        Anything layer 0 has mag 0 OR mag > 1/9e15 and < 9e15.
        Anything layer 1 or higher has abs(mag) >= 15.954 and < 9e15.
        We will assume in calculations that all Decimals are either erroneous or satisfy these criteria. (Otherwise: Garbage in, garbage out.)
        */
        if (this.sign === 0 || (this.mag === 0 && this.layer === 0))
        {
          this.sign = 0;
          this.mag = 0;
          this.layer = 0;
          return this;
        }
        
        if (this.layer === 0 && this.mag < 0)
        {
          //extract sign from negative mag at layer 0
          this.mag = -this.mag;
          this.sign = -this.sign;
        }
        
        //Handle shifting from layer 0 to negative layers.
        if (this.layer === 0 && this.mag < FIRST_NEG_LAYER)
        {
          this.layer += 1;
          this.mag = Math.log10(this.mag);
          return this;
        }
        
        var absmag = Math.abs(this.mag);
        var signmag = Math.sign(this.mag);
        
        if (absmag >= EXP_LIMIT)
        {
          this.layer += 1;
          this.mag = signmag*Math.log10(absmag);
          return this;
        }
        else
        {
          while (absmag < LAYER_DOWN && this.layer > 0)
          {
            this.layer -= 1;
            if (this.layer === 0)
            {
              this.mag = Math.pow(10, this.mag);
            }
            else
            {
              this.mag = signmag*Math.pow(10, absmag);
              absmag = Math.abs(this.mag);
              signmag = Math.sign(this.mag);
            }
          }
          if (this.layer === 0)
          {
            if (this.mag < 0)
            {
              //extract sign from negative mag at layer 0
              this.mag = -this.mag;
              this.sign = -this.sign;
            }
            else if (this.mag === 0)
            {
              //excessive rounding can give us all zeroes
              this.sign = 0;
            }
          }
        }
  
        return this;
      };
  
      Decimal.prototype.fromComponents = function (sign, layer, mag) {
        this.sign = sign;
        this.layer = layer;
        this.mag = mag;
  
        this.normalize();
        return this;
      };
  
      Decimal.prototype.fromComponents_noNormalize = function (sign, layer, mag) {
        this.sign = sign;
        this.layer = layer;
        this.mag = mag;
        return this;
      };
      
      Decimal.prototype.fromMantissaExponent = function (mantissa, exponent) {
        this.layer = 1;
        this.sign = Math.sign(mantissa);
        mantissa = Math.abs(mantissa);
        this.mag = exponent + Math.log10(mantissa);
  
        this.normalize();
        return this;
      };
  
  
      Decimal.prototype.fromMantissaExponent_noNormalize = function (mantissa, exponent) {
        //The idea of 'normalizing' a break_infinity.js style Decimal doesn't really apply. So just do the same thing.
        this.fromMantissaExponent(mantissa, exponent);
        return this;
      };
  
      Decimal.prototype.fromDecimal = function (value) {
        this.sign = value.sign;
        this.layer = value.layer;
        this.mag = value.mag;
        return this;
      };
  
      Decimal.prototype.fromNumber = function (value) {
        this.mag = Math.abs(value);
        this.sign = Math.sign(value);
        this.layer = 0;
        this.normalize();
        return this;
      };
  
      var IGNORE_COMMAS = true;
      var COMMAS_ARE_DECIMAL_POINTS = false;
      
      Decimal.prototype.fromString = function (value) {
        if (IGNORE_COMMAS) { value = value.replace(",", ""); }
        else if (COMMAS_ARE_DECIMAL_POINTS) { value = value.replace(",", "."); }
      
        //Handle x^^^y format.
        var pentationparts = value.split("^^^");
        if (pentationparts.length === 2)
        {
          var base = parseFloat(pentationparts[0]);
          var height = parseFloat(pentationparts[1]);
          var payload = 1;
          var heightparts = pentationparts[1].split(";");
          if (heightparts.length === 2)
          {
            var payload = parseFloat(heightparts[1]);
            if (!isFinite(payload)) { payload = 1; }
          }
          if (isFinite(base) && isFinite(height))
          {
            var result = Decimal.pentate(base, height, payload);
            this.sign = result.sign;
            this.layer = result.layer;
            this.mag = result.mag;
            return this;
          }
        }
      
        //Handle x^^y format.
        var tetrationparts = value.split("^^");
        if (tetrationparts.length === 2)
        {
          var base = parseFloat(tetrationparts[0]);
          var height = parseFloat(tetrationparts[1]);
          var heightparts = tetrationparts[1].split(";");
          if (heightparts.length === 2)
          {
            var payload = parseFloat(heightparts[1]);
            if (!isFinite(payload)) { payload = 1; }
          }
          if (isFinite(base) && isFinite(height))
          {
            var result = Decimal.tetrate(base, height, payload);
            this.sign = result.sign;
            this.layer = result.layer;
            this.mag = result.mag;
            return this;
          }
        }
        
        //Handle x^y format.
        var powparts = value.split("^");
        if (powparts.length === 2)
        {
          var base = parseFloat(powparts[0]);
          var exponent = parseFloat(powparts[1]);
          if (isFinite(base) && isFinite(exponent))
          {
            var result = Decimal.pow(base, exponent);
            this.sign = result.sign;
            this.layer = result.layer;
            this.mag = result.mag;
            return this;
          }
        }
        
        //Handle various cases involving it being a Big Number.
        value = value.trim().toLowerCase();
        
        //handle X PT Y format.
        var ptparts = value.split("pt");
        if (ptparts.length === 2)
        {
          base = 10;
          height = parseFloat(ptparts[0]);
          ptparts[1] = ptparts[1].replace("(", "");
          ptparts[1] = ptparts[1].replace(")", "");
          var payload = parseFloat(ptparts[1]);
          if (!isFinite(payload)) { payload = 1; }
          if (isFinite(base) && isFinite(height))
          {
            var result = Decimal.tetrate(base, height, payload);
            this.sign = result.sign;
            this.layer = result.layer;
            this.mag = result.mag;
            return this;
          }
        }
        
        //handle XpY format (it's the same thing just with p).
        var ptparts = value.split("p");
        if (ptparts.length === 2)
        {
          base = 10;
          height = parseFloat(ptparts[0]);
          ptparts[1] = ptparts[1].replace("(", "");
          ptparts[1] = ptparts[1].replace(")", "");
          var payload = parseFloat(ptparts[1]);
          if (!isFinite(payload)) { payload = 1; }
          if (isFinite(base) && isFinite(height))
          {
            var result = Decimal.tetrate(base, height, payload);
            this.sign = result.sign;
            this.layer = result.layer;
            this.mag = result.mag;
            return this;
          }
        }
  
        var parts = value.split("e");
        var ecount = parts.length-1;
      
        //Handle numbers that are exactly floats (0 or 1 es).
        if (ecount === 0)
        {
          var numberAttempt = parseFloat(value);
          if (isFinite(numberAttempt))
          {
            return this.fromNumber(numberAttempt);
          }
        }
        else if (ecount === 1)
        {
          //Very small numbers ("2e-3000" and so on) may look like valid floats but round to 0.
          var numberAttempt = parseFloat(value);
          if (isFinite(numberAttempt) && numberAttempt !== 0)
          {
            return this.fromNumber(numberAttempt);
          }
        }
        
        //Handle new (e^N)X format.
        var newparts = value.split("e^");
        if (newparts.length === 2)
        {
          this.sign = 1;
          if (newparts[0].charAt(0) == "-")
          {
            this.sign = -1;
          }
          var layerstring = "";
          for (var i = 0; i < newparts[1].length; ++i)
          {
            var chrcode = newparts[1].charCodeAt(i);
            if ((chrcode >= 43 && chrcode <= 57) || chrcode === 101) //is "0" to "9" or "+" or "-" or "." or "e" (or "," or "/")
            {
              layerstring += newparts[1].charAt(i);
            }
            else //we found the end of the layer count
            {
              this.layer = parseFloat(layerstring);
              this.mag = parseFloat(newparts[1].substr(i+1));
              this.normalize();
              return this;
            }
          }
        }
        
        if (ecount < 1) { this.sign = 0; this.layer = 0; this.mag = 0; return this; }
        var mantissa = parseFloat(parts[0]);
        if (mantissa === 0) { this.sign = 0; this.layer = 0; this.mag = 0; return this; }
        var exponent = parseFloat(parts[parts.length-1]);
        //handle numbers like AeBeC and AeeeeBeC
        if (ecount >= 2)
        {
          var me = parseFloat(parts[parts.length-2]);
          if (isFinite(me))
          {
            exponent *= Math.sign(me);
            exponent += f_maglog10(me);
          }
        }
        
        //Handle numbers written like eee... (N es) X
        if (!isFinite(mantissa))
        {
          this.sign = (parts[0] === "-") ? -1 : 1;
          this.layer = ecount;
          this.mag = exponent;
        }
        //Handle numbers written like XeY
        else if (ecount === 1)
        {
          this.sign = Math.sign(mantissa);
          this.layer = 1;
          //Example: 2e10 is equal to 10^log10(2e10) which is equal to 10^(10+log10(2))
          this.mag = exponent + Math.log10(Math.abs(mantissa));
        }
        //Handle numbers written like Xeee... (N es) Y
        else
        {
          this.sign = Math.sign(mantissa);
          this.layer = ecount;
          if (ecount === 2)
          {
            var result = Decimal.mul(FC(1, 2, exponent), D(mantissa));
            this.sign = result.sign;
            this.layer = result.layer;
            this.mag = result.mag;
            return this;
          }
          else
          {
            //at eee and above, mantissa is too small to be recognizable!
            this.mag = exponent;
          }
        }
        
        this.normalize();
        return this;
      };
  
      Decimal.prototype.fromValue = function (value) {
        if (value instanceof Decimal) {
          return this.fromDecimal(value);
        }
  
        if (typeof value === "number") {
          return this.fromNumber(value);
        }
  
        if (typeof value === "string") {
          return this.fromString(value);
        }
  
        this.sign = 0;
        this.layer = 0;
        this.mag = 0;
        return this;
      };
  
      Decimal.prototype.toNumber = function () {
        if (!Number.isFinite(this.layer)) { return Number.NaN; }
        if (this.layer === 0)
        {
          return this.sign*this.mag;
        }
        else if (this.layer === 1)
        {
          return this.sign*Math.pow(10, this.mag);
        }
        else //overflow for any normalized Decimal
        {
          return this.mag > 0 ? (this.sign > 0 ? Number.POSITIVE_INFINITY : Number.NEGATIVE_INFINITY) : 0;
        }
      };
      
      Decimal.prototype.mantissaWithDecimalPlaces = function (places) {
        // https://stackoverflow.com/a/37425022
        if (isNaN(this.m)) {
          return Number.NaN;
        }
  
        if (this.m === 0) {
          return 0;
        }
  
        return decimalPlaces(this.m, places);
      };
      
      Decimal.prototype.magnitudeWithDecimalPlaces = function (places) {
        // https://stackoverflow.com/a/37425022
        if (isNaN(this.mag)) {
          return Number.NaN;
        }
  
        if (this.mag === 0) {
          return 0;
        }
  
        return decimalPlaces(this.mag, places);
      };
      
      Decimal.prototype.toString = function () {
        if (this.layer === 0)
        {
          if ((this.mag < 1e21 && this.mag > 1e-7) || this.mag === 0)
          {
            return (this.sign*this.mag).toString();
          }
          return this.m + "e" + this.e;
        }
        else if (this.layer === 1)
        {
          return this.m + "e" + this.e;
        }
        else
        {
          //layer 2+
          if (this.layer <= MAX_ES_IN_A_ROW)
          {
            return (this.sign === -1 ? "-" : "") + "e".repeat(this.layer) + this.mag;
          }
          else
          {
            return (this.sign === -1 ? "-" : "") + "(e^" + this.layer + ")" + this.mag;
          }
        }
      };
      
      Decimal.prototype.toExponential = function (places) {
        if (this.layer === 0)
        {
          return (this.sign*this.mag).toExponential(places);
        }
        return this.toStringWithDecimalPlaces(places);
      };
      
      Decimal.prototype.toFixed = function (places) {
        if (this.layer === 0)
        {
          return (this.sign*this.mag).toFixed(places);
        }
        return this.toStringWithDecimalPlaces(places);
      };
      
      Decimal.prototype.toPrecision = function (places) {
        if (this.e <= -7) {
          return this.toExponential(places - 1);
        }
  
        if (places > this.e) {
          return this.toFixed(places - this.exponent - 1);
        }
  
        return this.toExponential(places - 1);
      };
      
      Decimal.prototype.valueOf = function () {
        return this.toString();
      };
  
      Decimal.prototype.toJSON = function () {
        return this.toString();
      };
      
      Decimal.prototype.toStringWithDecimalPlaces = function (places) {
        if (this.layer === 0)
        {
          if ((this.mag < 1e21 && this.mag > 1e-7) || this.mag === 0)
          {
            return (this.sign*this.mag).toFixed(places);
          }
          return decimalPlaces(this.m, places) + "e" + decimalPlaces(this.e, places);
        }
        else if (this.layer === 1)
        {
          return decimalPlaces(this.m, places) + "e" + decimalPlaces(this.e, places);
        }
        else
        {
          //layer 2+
          if (this.layer <= MAX_ES_IN_A_ROW)
          {
            return (this.sign === -1 ? "-" : "") + "e".repeat(this.layer) + decimalPlaces(this.mag, places);
          }
          else
          {
            return (this.sign === -1 ? "-" : "") + "(e^" + this.layer + ")" + decimalPlaces(this.mag, places);
          }
        }
      };
      
      Decimal.prototype.abs = function () {
        return FC_NN(this.sign === 0 ? 0 : 1, this.layer, this.mag);
      };
  
      Decimal.prototype.neg = function () {
        return FC_NN(-this.sign, this.layer, this.mag);
      };
  
      Decimal.prototype.negate = function () {
        return this.neg();
      };
  
      Decimal.prototype.negated = function () {
        return this.neg();
      };
  
      Decimal.prototype.sign = function () {
        return this.sign;
      };
  
      Decimal.prototype.sgn = function () {
        return this.sign;
      };
      
      Decimal.prototype.round = function () {
        if (this.mag < 0)
        {
          return Decimal.dZero;
        }
        if (this.layer === 0)
        {
          return FC(this.sign, 0, Math.round(this.mag));
        }
        return this;
      };
  
      Decimal.prototype.floor = function () {
        if (this.mag < 0)
        {
          return Decimal.dZero;
        }
        if (this.layer === 0)
        {
          return FC(this.sign, 0, Math.floor(this.mag));
        }
        return this;
      };
  
      Decimal.prototype.ceil = function () {
        if (this.mag < 0)
        {
          return Decimal.dZero;
        }
        if (this.layer === 0)
        {
          return FC(this.sign, 0, Math.ceil(this.mag));
        }
        return this;
      };
  
      Decimal.prototype.trunc = function () {
        if (this.mag < 0)
        {
          return Decimal.dZero;
        }
        if (this.layer === 0)
        {
          return FC(this.sign, 0, Math.trunc(this.mag));
        }
        return this;
      };
  
      Decimal.prototype.add = function (value) {
        var decimal = D(value);
        
        //inf/nan check
        if (!Number.isFinite(this.layer)) { return this; }
        if (!Number.isFinite(decimal.layer)) { return decimal; }
        
        //Special case - if one of the numbers is 0, return the other number.
        if (this.sign === 0) { return decimal; }
        if (decimal.sign === 0) { return this; }
        
        //Special case - Adding a number to its negation produces 0, no matter how large.
        if (this.sign === -(decimal.sign) && this.layer === decimal.layer && this.mag === decimal.mag) { return FC_NN(0, 0, 0); }
        
        var a;
        var b;
        
        //Special case: If one of the numbers is layer 2 or higher, just take the bigger number.
        if ((this.layer >= 2 || decimal.layer >= 2)) { return this.maxabs(decimal); }
        
        if (Decimal.cmpabs(this, decimal) > 0)
        {
          a = this;
          b = decimal;
        }
        else
        {
          a = decimal;
          b = this;
        }
        
        if (a.layer === 0 && b.layer === 0) { return D(a.sign*a.mag + b.sign*b.mag); }
        
        var layera = a.layer*Math.sign(a.mag);
        var layerb = b.layer*Math.sign(b.mag);
        
        //If one of the numbers is 2+ layers higher than the other, just take the bigger number.
        if (layera - layerb >= 2) { return a; }
        
        if (layera === 0 && layerb === -1)
        {
          if (Math.abs(b.mag-Math.log10(a.mag)) > MAX_SIGNIFICANT_DIGITS)
          {
            return a;
          }
          else
          {
            var magdiff = Math.pow(10, Math.log10(a.mag)-b.mag);
            var mantissa = (b.sign)+(a.sign*magdiff);
            return FC(Math.sign(mantissa), 1, b.mag+Math.log10(Math.abs(mantissa)));
          }
        }
        
        if (layera === 1 && layerb === 0)
        {
          if (Math.abs(a.mag-Math.log10(b.mag)) > MAX_SIGNIFICANT_DIGITS)
          {
            return a;
          }
          else
          {
            var magdiff = Math.pow(10, a.mag-Math.log10(b.mag));
            var mantissa = (b.sign)+(a.sign*magdiff);
            return FC(Math.sign(mantissa), 1, Math.log10(b.mag)+Math.log10(Math.abs(mantissa)));
          }
        }
        
        if (Math.abs(a.mag-b.mag) > MAX_SIGNIFICANT_DIGITS)
        {
          return a;
        }
        else
        {
          var magdiff = Math.pow(10, a.mag-b.mag);
          var mantissa = (b.sign)+(a.sign*magdiff);
          return FC(Math.sign(mantissa), 1, b.mag+Math.log10(Math.abs(mantissa)));
        }
        
        throw Error("Bad arguments to add: " + this + ", " + value);
      };
  
      Decimal.prototype.plus = function (value) {
        return this.add(value);
      };
  
      Decimal.prototype.sub = function (value) {
        return this.add(D(value).neg());
      };
  
      Decimal.prototype.subtract = function (value) {
        return this.sub(value);
      };
  
      Decimal.prototype.minus = function (value) {
        return this.sub(value);
      };
  
      Decimal.prototype.mul = function (value) {
        var decimal = D(value);
        
        //inf/nan check
        if (!Number.isFinite(this.layer)) { return this; }
        if (!Number.isFinite(decimal.layer)) { return decimal; }
        
        //Special case - if one of the numbers is 0, return 0.
        if (this.sign === 0 || decimal.sign === 0) { return FC_NN(0, 0, 0); }
        
        //Special case - Multiplying a number by its own reciprocal yields +/- 1, no matter how large.
        if (this.layer === decimal.layer && this.mag === -decimal.mag) { return FC_NN(this.sign*decimal.sign, 0, 1); }
              
        var a;
        var b;
        
        //Which number is bigger in terms of its multiplicative distance from 1?
        if ((this.layer > decimal.layer) || (this.layer == decimal.layer && Math.abs(this.mag) > Math.abs(decimal.mag)))
        {
          a = this;
          b = decimal;
        }
        else
        {
          a = decimal;
          b = this;
        }
        
        if (a.layer === 0 && b.layer === 0) { return D(a.sign*b.sign*a.mag*b.mag); }
        
        //Special case: If one of the numbers is layer 3 or higher or one of the numbers is 2+ layers bigger than the other, just take the bigger number.
        if (a.layer >= 3 || (a.layer - b.layer >= 2)) { return FC(a.sign*b.sign, a.layer, a.mag); }
  
        if (a.layer === 1 && b.layer === 0)
        { 
          return FC(a.sign*b.sign, 1, a.mag+Math.log10(b.mag));
        }
        
        if (a.layer === 1 && b.layer === 1)
        {
          return FC(a.sign*b.sign, 1, a.mag+b.mag);
        }
        
        if (a.layer === 2 && b.layer === 1)
        {
          var newmag = FC(Math.sign(a.mag), a.layer-1, Math.abs(a.mag)).add(FC(Math.sign(b.mag), b.layer-1, Math.abs(b.mag)));
          return FC(a.sign*b.sign, newmag.layer+1, newmag.sign*newmag.mag);
        }
        
        if (a.layer === 2 && b.layer === 2)
        {
          var newmag = FC(Math.sign(a.mag), a.layer-1, Math.abs(a.mag)).add(FC(Math.sign(b.mag), b.layer-1, Math.abs(b.mag)));
          return FC(a.sign*b.sign, newmag.layer+1, newmag.sign*newmag.mag);
        }
        
        throw Error("Bad arguments to mul: " + this + ", " + value);
      };
  
      Decimal.prototype.multiply = function (value) {
        return this.mul(value);
      };
  
      Decimal.prototype.times = function (value) {
        return this.mul(value);
      };
  
      Decimal.prototype.div = function (value) {
        var decimal = D(value);
        return this.mul(decimal.recip());
      };
  
      Decimal.prototype.divide = function (value) {
        return this.div(value);
      };
  
      Decimal.prototype.divideBy = function (value) {
        return this.div(value);
      };
  
      Decimal.prototype.dividedBy = function (value) {
        return this.div(value);
      };
  
      Decimal.prototype.recip = function () {
        if (this.mag === 0)
        {
          return Decimal.dNaN;
        }
        else if (this.layer === 0)
        {
          return FC(this.sign, 0, 1/this.mag);
        }
        else
        {
          return FC(this.sign, this.layer, -this.mag);
        }
      };
  
      Decimal.prototype.reciprocal = function () {
        return this.recip();
      };
  
      Decimal.prototype.reciprocate = function () {
        return this.recip();
      };
      
      /**
       * -1 for less than value, 0 for equals value, 1 for greater than value
       */
      Decimal.prototype.cmp = function (value) {
        var decimal = D(value);
        if (this.sign > decimal.sign) { return 1; }
        if (this.sign < decimal.sign) { return -1; }
        return this.sign*this.cmpabs(value);
      };
      
      Decimal.prototype.cmpabs = function (value) {
        var decimal = D(value);
        var layera = this.mag > 0 ? this.layer : -this.layer;
        var layerb = decimal.mag > 0 ? decimal.layer : -decimal.layer;
        if (layera > layerb) { return 1; }
        if (layera < layerb) { return -1; }
        if (this.mag > decimal.mag) { return 1; }
        if (this.mag < decimal.mag) { return -1; }
        return 0;
      };
  
      Decimal.prototype.compare = function (value) {
        return this.cmp(value);
      };
  
      Decimal.prototype.eq = function (value) {
        var decimal = D(value);
        return this.sign === decimal.sign && this.layer === decimal.layer && this.mag === decimal.mag;
      };
  
      Decimal.prototype.equals = function (value) {
        return this.eq(value);
      };
  
      Decimal.prototype.neq = function (value) {
        return !this.eq(value);
      };
  
      Decimal.prototype.notEquals = function (value) {
        return this.neq(value);
      };
  
      Decimal.prototype.lt = function (value) {
        var decimal = D(value);
        return this.cmp(value) === -1;
      };
  
      Decimal.prototype.lte = function (value) {
        return !this.gt(value);
      };
  
      Decimal.prototype.gt = function (value) {
        var decimal = D(value);
        return this.cmp(value) === 1;
      };
  
      Decimal.prototype.gte = function (value) {
        return !this.lt(value);
      };
  
      Decimal.prototype.max = function (value) {
        var decimal = D(value);
        return this.lt(decimal) ? decimal : this;
      };
  
      Decimal.prototype.min = function (value) {
        var decimal = D(value);
        return this.gt(decimal) ? decimal : this;
      };
      
      Decimal.prototype.maxabs = function (value) {
        var decimal = D(value);
        return this.cmpabs(decimal) < 0 ? decimal : this;
      };
  
      Decimal.prototype.minabs = function (value) {
        var decimal = D(value);
        return this.cmpabs(decimal) > 0 ? decimal : this;
      };
      
      Decimal.prototype.clamp = function(min, max) {
        return this.max(min).min(max);
      }
      
      Decimal.prototype.clampMin = function(min) {
        return this.max(min);
      }
      
      Decimal.prototype.clampMax = function(max) {
        return this.min(max);
      }
  
      Decimal.prototype.cmp_tolerance = function (value, tolerance) {
        var decimal = D(value);
        return this.eq_tolerance(decimal, tolerance) ? 0 : this.cmp(decimal);
      };
  
      Decimal.prototype.compare_tolerance = function (value, tolerance) {
        return this.cmp_tolerance(value, tolerance);
      };
      
      /**
       * Tolerance is a relative tolerance, multiplied by the greater of the magnitudes of the two arguments.
       * For example, if you put in 1e-9, then any number closer to the
       * larger number than (larger number)*1e-9 will be considered equal.
       */
      Decimal.prototype.eq_tolerance = function (value, tolerance) {
        var decimal = D(value); // https://stackoverflow.com/a/33024979
        if (tolerance == null) { tolerance = 1e-7; }
        //Numbers that are too far away are never close.
        if (this.sign !== decimal.sign) { return false; }
        if (Math.abs(this.layer - decimal.layer) > 1) { return false; }
        // return abs(a-b) <= tolerance * max(abs(a), abs(b))
        var magA = this.mag;
        var magB = decimal.mag;
        if (this.layer > decimal.layer) { magB = f_maglog10(magB); }
        if (this.layer < decimal.layer) { magA = f_maglog10(magA); }
        return Math.abs(magA-magB) <= tolerance*Math.max(Math.abs(magA), Math.abs(magB));
      };
  
      Decimal.prototype.equals_tolerance = function (value, tolerance) {
        return this.eq_tolerance(value, tolerance);
      };
  
      Decimal.prototype.neq_tolerance = function (value, tolerance) {
        return !this.eq_tolerance(value, tolerance);
      };
  
      Decimal.prototype.notEquals_tolerance = function (value, tolerance) {
        return this.neq_tolerance(value, tolerance);
      };
  
      Decimal.prototype.lt_tolerance = function (value, tolerance) {
        var decimal = D(value);
        return !this.eq_tolerance(decimal, tolerance) && this.lt(decimal);
      };
  
      Decimal.prototype.lte_tolerance = function (value, tolerance) {
        var decimal = D(value);
        return this.eq_tolerance(decimal, tolerance) || this.lt(decimal);
      };
  
      Decimal.prototype.gt_tolerance = function (value, tolerance) {
        var decimal = D(value);
        return !this.eq_tolerance(decimal, tolerance) && this.gt(decimal);
      };
  
      Decimal.prototype.gte_tolerance = function (value, tolerance) {
        var decimal = D(value);
        return this.eq_tolerance(decimal, tolerance) || this.gt(decimal);
      };
      
      Decimal.prototype.pLog10 = function() {
        if (this.lt(Decimal.dZero)) { return Decimal.dZero; }
        return this.log10();
      }
  
      Decimal.prototype.absLog10 = function () {
        if (this.sign === 0)
        {
          return Decimal.dNaN;
        }
        else if (this.layer > 0)
        {
          return FC(Math.sign(this.mag), this.layer-1, Math.abs(this.mag));
        }
        else
        {
          return FC(1, 0, Math.log10(this.mag));
        }
      };
      
      Decimal.prototype.log10 = function () {
        if (this.sign <= 0)
        {
          return Decimal.dNaN;
        }
        else if (this.layer > 0)
        {
          return FC(Math.sign(this.mag), this.layer-1, Math.abs(this.mag));
        }
        else
        {
          return FC(this.sign, 0, Math.log10(this.mag));
        }
      };
  
      Decimal.prototype.log = function (base) {
        base = D(base);
        if (this.sign <= 0)
        {
          return Decimal.dNaN;
        }
        if (base.sign <= 0)
        {
          return Decimal.dNaN;
        }
        if (base.sign === 1 && base.layer === 0 && base.mag === 1)
        {
          return Decimal.dNaN;
        }
        else if (this.layer === 0 && base.layer === 0)
        {
          return FC(this.sign, 0, Math.log(this.mag)/Math.log(base.mag));
        }
        
        return Decimal.div(this.log10(), base.log10());
      };
  
      Decimal.prototype.log2 = function () {
        if (this.sign <= 0)
        {
          return Decimal.dNaN;
        }
        else if (this.layer === 0)
        {
          return FC(this.sign, 0, Math.log2(this.mag));
        }
        else if (this.layer === 1)
        {
          return FC(Math.sign(this.mag), 0, Math.abs(this.mag)*3.321928094887362); //log2(10)
        }
        else if (this.layer === 2)
        {
          return FC(Math.sign(this.mag), 1, Math.abs(this.mag)+0.5213902276543247); //-log10(log10(2))
        }
        else
        {
          return FC(Math.sign(this.mag), this.layer-1, Math.abs(this.mag));
        }
      };
  
      Decimal.prototype.ln = function () {
        if (this.sign <= 0)
        {
          return Decimal.dNaN;
        }
        else if (this.layer === 0)
        {
          return FC(this.sign, 0, Math.log(this.mag));
        }
        else if (this.layer === 1)
        {
          return FC(Math.sign(this.mag), 0, Math.abs(this.mag)*2.302585092994046); //ln(10)
        }
        else if (this.layer === 2)
        {
          return FC(Math.sign(this.mag), 1, Math.abs(this.mag)+0.36221568869946325); //log10(log10(e))
        }
        else
        {
          return FC(Math.sign(this.mag), this.layer-1, Math.abs(this.mag));
        }
      };
  
      Decimal.prototype.logarithm = function (base) {
        return this.log(base);
      };
  
      Decimal.prototype.pow = function (value) {
        var decimal = D(value);
        var a = this;
        var b = decimal;
  
        //special case: if a is 0, then return 0 (UNLESS b is 0, then return 1)
        if (a.sign === 0) { return b.eq(0) ? FC_NN(1, 0, 1) : a; }
        //special case: if a is 1, then return 1
        if (a.sign === 1 && a.layer === 0 && a.mag === 1) { return a; }
        //special case: if b is 0, then return 1
        if (b.sign === 0) { return FC_NN(1, 0, 1); }
        //special case: if b is 1, then return a
        if (b.sign === 1 && b.layer === 0 && b.mag === 1) { return a; }
        
        var result = (a.absLog10().mul(b)).pow10();
  
        if (this.sign === -1 && Math.abs(b.toNumber() % 2) === 1) {
          return result.neg();
        }
  
        return result;
      };
      
      Decimal.prototype.pow10 = function() {
        /*
        There are four cases we need to consider:
        1) positive sign, positive mag (e15, ee15): +1 layer (e.g. 10^15 becomes e15, 10^e15 becomes ee15)
        2) negative sign, positive mag (-e15, -ee15): +1 layer but sign and mag sign are flipped (e.g. 10^-15 becomes e-15, 10^-e15 becomes ee-15)
        3) positive sign, negative mag (e-15, ee-15): layer 0 case would have been handled in the Math.pow check, so just return 1
        4) negative sign, negative mag (-e-15, -ee-15): layer 0 case would have been handled in the Math.pow check, so just return 1
        */
        
        if (!Number.isFinite(this.layer) || !Number.isFinite(this.mag)) { return Decimal.dNaN; }
        
        var a = this;
        
        //handle layer 0 case - if no precision is lost just use Math.pow, else promote one layer
        if (a.layer === 0)
        {
          var newmag = Math.pow(10, a.sign*a.mag);
          if (Number.isFinite(newmag) && Math.abs(newmag) > 0.1) { return FC(1, 0, newmag); }
          else
          {
            if (a.sign === 0) { return Decimal.dOne; }
            else { a = FC_NN(a.sign, a.layer+1, Math.log10(a.mag)); }
          }
        }
        
        //handle all 4 layer 1+ cases individually
        if (a.sign > 0 && a.mag > 0)
        {
          return FC(a.sign, a.layer+1, a.mag);
        }
        if (a.sign < 0 && a.mag > 0)
        {
          return FC(-a.sign, a.layer+1, -a.mag);
        }
        //both the negative mag cases are identical: one +/- rounding error
        return Decimal.dOne;
      }
  
      Decimal.prototype.pow_base = function (value) {
        return D(value).pow(this);
      };
      
      Decimal.prototype.root = function (value) {
        var decimal = D(value);
        return this.pow(decimal.recip());
      }
  
      Decimal.prototype.factorial = function () {
        if (this.mag < 0)
        {
          return this.toNumber().add(1).gamma();
        }
        else if (this.layer === 0)
        {
          return this.add(1).gamma();
        }
        else if (this.layer === 1)
        {
          return Decimal.exp(Decimal.mul(this, Decimal.ln(this).sub(1)));
        }
        else
        {
          return Decimal.exp(this);
        }
      };
      
      //from HyperCalc source code
      Decimal.prototype.gamma = function () {
        if (this.mag < 0)
        {
          return this.recip();
        }
        else if (this.layer === 0)
        {
          if (this.lt(FC_NN(1, 0, 24)))
          {
            return D(f_gamma(this.sign*this.mag));
          }
          
          var t = this.mag - 1;
          var l = 0.9189385332046727; //0.5*Math.log(2*Math.PI)
          l = (l+((t+0.5)*Math.log(t)));
          l = l-t;
          var n2 = t*t;
          var np = t;
          var lm = 12*np;
          var adj = 1/lm;
          var l2 = l+adj;
          if (l2 === l)
          {
            return Decimal.exp(l);
          }
          
          l = l2;
          np = np*n2;
          lm = 360*np;
          adj = 1/lm;
          l2 = l-adj;
          if (l2 === l)
          {
            return Decimal.exp(l);
          }
          
          l = l2;
          np = np*n2;
          lm = 1260*np;
          var lt = 1/lm;
          l = l+lt;
          np = np*n2;
          lm = 1680*np;
          lt = 1/lm;
          l = l-lt;
          return Decimal.exp(l);
        }
        else if (this.layer === 1)
        {
          return Decimal.exp(Decimal.mul(this, Decimal.ln(this).sub(1)));
        }
        else
        {
          return Decimal.exp(this);
        }
      };
      
      Decimal.prototype.lngamma = function () {
        return this.gamma().ln();
      }
  
      Decimal.prototype.exp = function () {
        if (this.mag < 0) { return Decimal.dOne; }
        if (this.layer === 0 && this.mag <= 709.7) { return D(Math.exp(this.sign*this.mag)); }
        else if (this.layer === 0) { return FC(1, 1, this.sign*Math.log10(Math.E)*this.mag); }
        else if (this.layer === 1) { return FC(1, 2, this.sign*(Math.log10(0.4342944819032518)+this.mag)); }
        else { return FC(1, this.layer+1, this.sign*this.mag); }
      };
  
      Decimal.prototype.sqr = function () {
        return this.pow(2);
      };
  
      Decimal.prototype.sqrt = function () {
        if (this.layer === 0) { return D(Math.sqrt(this.sign*this.mag)); }
        else if (this.layer === 1) { return FC(1, 2, Math.log10(this.mag)-0.3010299956639812); }
        else
        {
          var result = Decimal.div(FC_NN(this.sign, this.layer-1, this.mag), FC_NN(1, 0, 2));
          result.layer += 1;
          result.normalize();
          return result;
        }
      };
  
      Decimal.prototype.cube = function () {
        return this.pow(3);
      };
  
      Decimal.prototype.cbrt = function () {
        return this.pow(1/3);
      };
      
      //Tetration/tetrate: The result of exponentiating 'this' to 'this' 'height' times in a row.  https://en.wikipedia.org/wiki/Tetration
      //If payload != 1, then this is 'iterated exponentiation', the result of exping (payload) to base (this) (height) times. https://andydude.github.io/tetration/archives/tetration2/ident.html
      //Works with negative and positive real heights.
      Decimal.prototype.tetrate = function(height = 2, payload = FC_NN(1, 0, 1)) {
        if (height === Number.POSITIVE_INFINITY)
        {
          //Formula for infinite height power tower.
          var negln = Decimal.ln(this).neg();
          return negln.lambertw().div(negln);
        }
        
        if (height < 0)
        {
          return Decimal.iteratedlog(payload, this, -height);
        }
        
        payload = D(payload);
        var oldheight = height;
        height = Math.trunc(height);
        var fracheight = oldheight-height;
       
        if (fracheight !== 0)
        {
          if (payload.eq(Decimal.dOne))
          {
            ++height;
            payload = new Decimal(fracheight);
          }
          else
          {
            if (this.eq(10))
            {
              payload = payload.layeradd10(fracheight);
            }
            else
            {
              payload = payload.layeradd(fracheight, this);
            }
          }
        }
        
        for (var i = 0; i < height; ++i)
        {
          payload = this.pow(payload);
          //bail if we're NaN
          if (!isFinite(payload.layer) || !isFinite(payload.mag)) { return payload; }
          //shortcut 
          if (payload.layer - this.layer > 3) { return FC_NN(payload.sign, payload.layer + (height - i - 1), payload.mag); }
          //give up after 100 iterations if nothing is happening
          if (i > 100) { return payload; }
        }
        return payload;
      }
      
      //iteratedexp/iterated exponentiation: - all cases handled in tetrate, so just call it
      Decimal.prototype.iteratedexp = function(height = 2, payload = FC_NN(1, 0, 1)) {
        return this.tetrate(height, payload);
      }
      
      //iterated log/repeated log: The result of applying log(base) 'times' times in a row. Approximately equal to subtracting (times) from the number's slog representation. Equivalent to tetrating to a negative height.
      //Works with negative and positive real heights.
      Decimal.prototype.iteratedlog = function(base = 10, times = 1) {      
        if (times < 0)
        {
          return Decimal.tetrate(base, -times, this);
        }
        
        base = D(base);
        var result = D(this);
        var fulltimes = times;
        times = Math.trunc(times);
        var fraction = fulltimes - times;
        if (result.layer - base.layer > 3)
        {
          var layerloss = Math.min(times, (result.layer - base.layer - 3));
          times -= layerloss;
          result.layer -= layerloss;
        }
        
        for (var i = 0; i < times; ++i)
        {
          result = result.log(base);
          //bail if we're NaN
          if (!isFinite(result.layer) || !isFinite(result.mag)) { return result; }
          //give up after 100 iterations if nothing is happening
          if (i > 100) { return result; }
        }
        
        //handle fractional part
        if (fraction > 0 && fraction < 1)
        {
          if (base.eq(10))
          {
            result = result.layeradd10(-fraction);
          }
          else
          {
            result = result.layeradd(-fraction, base);
          }
        }
        
        return result;
      }
      
      //Super-logarithm, one of tetration's inverses, tells you what size power tower you'd have to tetrate base to to get number. By definition, will never be higher than 1.8e308 in break_eternity.js, since a power tower 1.8e308 numbers tall is the largest representable number.
      // https://en.wikipedia.org/wiki/Super-logarithm
      Decimal.prototype.slog = function(base = 10) {
        if (this.mag < 0) { return Decimal.dNegOne; }
        
        base = D(base);
        
        var result = 0;
        var copy = D(this);
        if (copy.layer - base.layer > 3)
        {
          var layerloss = (copy.layer - base.layer - 3);
          result += layerloss;
          copy.layer -= layerloss;
        }
        
        for (var i = 0; i < 100; ++i)
        {
          if (copy.lt(Decimal.dZero))
          {
            copy = Decimal.pow(base, copy);
            result -= 1;
          }
          else if (copy.lte(Decimal.dOne))
          {
            return D(result + copy.toNumber() - 1); //<-- THIS IS THE CRITICAL FUNCTION
            //^ Also have to change tetrate payload handling and layeradd10 if this is changed!
          }
          else
          {
            result += 1;
            copy = Decimal.log(copy, base);
          }
        }
        return D(result);
      }
      
      //Approximations taken from the excellent paper https://web.archive.org/web/20090201164836/http://tetration.itgo.com/paper.html !
      //Not using for now unless I can figure out how to use it in all the related functions.
      /*var slog_criticalfunction_1 = function(x, z) {
        z = z.toNumber();
        return -1 + z;
      }
      
      var slog_criticalfunction_2 = function(x, z) {
        z = z.toNumber();
        var lnx = x.ln();
        if (lnx.layer === 0)
        {
          lnx = lnx.toNumber();
          return -1 + z*2*lnx/(1+lnx) - z*z*(1-lnx)/(1+lnx);
        }
        else
        {
          var term1 = lnx.mul(z*2).div(lnx.add(1));
          var term2 = Decimal.sub(1, lnx).mul(z*z).div(lnx.add(1));
          Decimal.dNegOne.add(Decimal.sub(term1, term2));
        }
      }
      
      var slog_criticalfunction_3 = function(x, z) {
        z = z.toNumber();
        var lnx = x.ln();
        var lnx2 = lnx.sqr();
        var lnx3 = lnx.cube();
        if (lnx.layer === 0 && lnx2.layer === 0 && lnx3.layer === 0)
        {
          lnx = lnx.toNumber();
          lnx2 = lnx2.toNumber();
          lnx3 = lnx3.toNumber();
          
          var term1 = 6*z*(lnx+lnx3);
          var term2 = 3*z*z*(3*lnx2-2*lnx3);
          var term3 = 2*z*z*z*(1-lnx-2*lnx2+lnx3);
          var top = term1+term2+term3;
          var bottom = 2+4*lnx+5*lnx2+2*lnx3;
          
          return -1 + top/bottom;
        }
        else
        {
          var term1 = (lnx.add(lnx3)).mul(6*z);
          var term2 = (lnx2.mul(3).sub(lnx3.mul(2))).mul(3*z*z);
          var term3 = (Decimal.dOne.sub(lnx).sub(lnx2.mul(2)).add(lnx3)).mul(2*z*z*z);
          var top = term1.add(term2).add(term3);
          var bottom = new Decimal(2).add(lnx.mul(4)).add(lnx2.mul(5)).add(lnx3.mul(2));
          
          return Decimal.dNegOne.add(top.div(bottom));
        }
      }*/
      
      //Function for adding/removing layers from a Decimal, even fractional layers (e.g. its slog10 representation).
      //Everything continues to use the linear approximation ATM.
      Decimal.prototype.layeradd10 = function(diff) {
        diff = Decimal.fromValue_noAlloc(diff).toNumber();
        var result = D(this);
        if (diff >= 1)
        {
          var layeradd = Math.trunc(diff);
          diff -= layeradd;
          result.layer += layeradd;
        }
        if (diff <= -1)
        {
          var layeradd = Math.trunc(diff);
          diff -= layeradd;
          result.layer += layeradd;
          if (result.layer < 0)
          {
            for (var i = 0; i < 100; ++i)
            {
              result.layer++;
              result.mag = Math.log10(result.mag);
              if (!isFinite(result.mag)) { return result; }
              if (result.layer >= 0) { break; }
            }
          }
        }
        
        //layeradd10: like adding 'diff' to the number's slog(base) representation. Very similar to tetrate base 10 and iterated log base 10. Also equivalent to adding a fractional amount to the number's layer in its break_eternity.js representation.
        if (diff > 0)
        {
          var subtractlayerslater = 0;
          //Ironically, this edge case would be unnecessary if we had 'negative layers'.
          while (Number.isFinite(result.mag) && result.mag < 10)
          {
            result.mag = Math.pow(10, result.mag);
            ++subtractlayerslater;
          }
          
          //A^(10^B) === C, solve for B
          //B === log10(logA(C))
          
          if (result.mag > 1e10)
          {
            result.mag = Math.log10(result.mag);
            result.layer++;
          }
          
          //Note that every integer slog10 value, the formula changes, so if we're near such a number, we have to spend exactly enough layerdiff to hit it, and then use the new formula.
          var diffToNextSlog = Math.log10(Math.log(1e10)/Math.log(result.mag), 10);
          if (diffToNextSlog < diff)
          {
            result.mag = Math.log10(1e10);
            result.layer++;
            diff -= diffToNextSlog;
          }
          
          result.mag = Math.pow(result.mag, Math.pow(10, diff));
          
          while (subtractlayerslater > 0)
          {
            result.mag = Math.log10(result.mag);
            --subtractlayerslater;
          }
        }
        else if (diff < 0)
        {
          var subtractlayerslater = 0;
          
          while (Number.isFinite(result.mag) && result.mag < 10)
          {
            result.mag = Math.pow(10, result.mag);
            ++subtractlayerslater;
          }
          
          if (result.mag > 1e10)
          {
            result.mag = Math.log10(result.mag);
            result.layer++;
          }
          
          var diffToNextSlog = Math.log10(1/Math.log10(result.mag));
          if (diffToNextSlog > diff)
          {
            result.mag = 1e10;
            result.layer--;
            diff -= diffToNextSlog;
          }
          
          result.mag = Math.pow(result.mag, Math.pow(10, diff));
          
          while (subtractlayerslater > 0)
          {
            result.mag = Math.log10(result.mag);
            --subtractlayerslater;
          }
        }
        
        while (result.layer < 0)
        {
          result.layer++;
          result.mag = Math.log10(result.mag);
        }
        result.normalize();
        return result;
      }
      
      //layeradd: like adding 'diff' to the number's slog(base) representation. Very similar to tetrate base 'base' and iterated log base 'base'.
      Decimal.prototype.layeradd = function(diff, base) {
        var slogthis = this.slog(base).toNumber();
        var slogdest = slogthis+diff;
        if (slogdest >= 0)
        {
          return Decimal.tetrate(base, slogdest);
        }
        else if (!Number.isFinite(slogdest))
        {
          return Decimal.dNaN;
        }
        else if (slogdest >= -1)
        {
          return Decimal.log(Decimal.tetrate(base, slogdest+1), base);
        }
        else
        {
          Decimal.log(Decimal.log(Decimal.tetrate(base, slogdest+2), base), base);
        }
      }
      
      //The Lambert W function, also called the omega function or product logarithm, is the solution W(x) === x*e^x.
      // https://en.wikipedia.org/wiki/Lambert_W_function
      //Some special values, for testing: https://en.wikipedia.org/wiki/Lambert_W_function#Special_values
      Decimal.prototype.lambertw = function() {
        if (this.lt(-0.3678794411710499))
        {
          throw Error("lambertw is unimplemented for results less than -1, sorry!");
        }
        else if (this.mag < 0)
        {
          return D(f_lambertw(this.toNumber()));
        }
        else if (this.layer === 0)
        {
          return D(f_lambertw(this.sign*this.mag));
        }
        else if (this.layer === 1)
        {
          return d_lambertw(this);
        }
        else if (this.layer === 2)
        {
          return d_lambertw(this);
        }
        if (this.layer >= 3)
        {
          return FC_NN(this.sign, this.layer-1, this.mag);
        }
      }
    
      //from https://github.com/scipy/scipy/blob/8dba340293fe20e62e173bdf2c10ae208286692f/scipy/special/lambertw.pxd
      // The evaluation can become inaccurate very close to the branch point
      // at ``-1/e``. In some corner cases, `lambertw` might currently
      // fail to converge, or can end up on the wrong branch.
      var d_lambertw = function(z, tol = 1e-10) {
      var w;
      var ew, wew, wewz, wn;
      
      if (!Number.isFinite(z.mag)) { return z; }
      if (z === 0)
      {
        return z;
      }
      if (z === 1)
      {
        //Split out this case because the asymptotic series blows up
        return OMEGA;
      }
      
      var absz = Decimal.abs(z);
      //Get an initial guess for Halley's method
      w = Decimal.ln(z);
      
      //Halley's method; see 5.9 in [1]
      
      for (var i = 0; i < 100; ++i)
      {
        ew = Decimal.exp(-w);
        wewz = w.sub(z.mul(ew));
        wn = w.sub(wewz.div(w.add(1).sub((w.add(2)).mul(wewz).div((Decimal.mul(2, w).add(2))))));
        if (Decimal.abs(wn.sub(w)).lt(Decimal.abs(wn).mul(tol)))
        {
          return wn;
        }
        else
        {
          w = wn;
        }
      }
      
      throw Error("Iteration failed to converge: " + z);
      //return Decimal.dNaN;
      }
      
      //The super square-root function - what number, tetrated to height 2, equals this?
      //Other sroots are possible to calculate probably through guess and check methods, this one is easy though.
      // https://en.wikipedia.org/wiki/Tetration#Super-root
      Decimal.prototype.ssqrt = function() {
        if (this.sign == 1 && this.layer >= 3)
        {
            return FC_NN(this.sign, this.layer-1, this.mag)
        }
        var lnx = this.ln();
        return lnx.div(lnx.lambertw());
      }
  /*
  
  Unit tests for tetrate/iteratedexp/iteratedlog/layeradd10/layeradd/slog:
  
  for (var i = 0; i < 1000; ++i)
  {
      var first = Math.random()*100;
      var both = Math.random()*100;
      var expected = first+both+1;
      var result = new Decimal(10).layeradd10(first).layeradd10(both).slog();
      if (Number.isFinite(result.mag) && !Decimal.eq_tolerance(expected, result))
      {
          console.log(first + ", " + both);
      }
  }
  
  for (var i = 0; i < 1000; ++i)
  {
      var first = Math.random()*100;
      var both = Math.random()*100;
      first += both;
      var expected = first-both+1;
      var result = new Decimal(10).layeradd10(first).layeradd10(-both).slog();
      if (Number.isFinite(result.mag) && !Decimal.eq_tolerance(expected, result))
      {
          console.log(first + ", " + both);
      }
  }
  
  for (var i = 0; i < 1000; ++i)
  {
      var first = Math.random()*100;
      var both = Math.random()*100;
      var base = Math.random()*8+2;
      var expected = first+both+1;
      var result = new Decimal(base).layeradd(first, base).layeradd(both, base).slog(base);
      if (Number.isFinite(result.mag) && !Decimal.eq_tolerance(expected, result))
      {
          console.log(first + ", " + both);
      }
  }
  
  for (var i = 0; i < 1000; ++i)
  {
      var first = Math.random()*100;
      var both = Math.random()*100;
      var base = Math.random()*8+2;
      first += both;
      var expected = first-both+1;
      var result = new Decimal(base).layeradd(first, base).layeradd(-both, base).slog(base);
      if (Number.isFinite(result.mag) && !Decimal.eq_tolerance(expected, result))
      {
          console.log(first + ", " + both);
      }
  }
  
  for (var i = 0; i < 1000; ++i)
  {
      var first = Math.round((Math.random()*30))/10;
      var both = Math.round((Math.random()*30))/10;
      var tetrateonly = Decimal.tetrate(10, first);
      var tetrateandlog = Decimal.tetrate(10, first+both).iteratedlog(10, both);
      if (!Decimal.eq_tolerance(tetrateonly, tetrateandlog))
      {
          console.log(first + ", " + both);
      }
  }
  
  for (var i = 0; i < 1000; ++i)
  {
      var first = Math.round((Math.random()*30))/10;
      var both = Math.round((Math.random()*30))/10;
    var base = Math.random()*8+2;
      var tetrateonly = Decimal.tetrate(base, first);
      var tetrateandlog = Decimal.tetrate(base, first+both).iteratedlog(base, both);
      if (!Decimal.eq_tolerance(tetrateonly, tetrateandlog))
      {
          console.log(first + ", " + both);
      }
  }
  
  for (var i = 0; i < 1000; ++i)
  {
      var first = Math.round((Math.random()*30))/10;
      var both = Math.round((Math.random()*30))/10;
    var base = Math.random()*8+2;
      var tetrateonly = Decimal.tetrate(base, first, base);
      var tetrateandlog = Decimal.tetrate(base, first+both, base).iteratedlog(base, both);
      if (!Decimal.eq_tolerance(tetrateonly, tetrateandlog))
      {
          console.log(first + ", " + both);
      }
  }
  
  for (var i = 0; i < 1000; ++i)
  {
      var xex = new Decimal(-0.3678794411710499+Math.random()*100);
      var x = Decimal.lambertw(xex);
      if (!Decimal.eq_tolerance(xex, x.mul(Decimal.exp(x))))
      {
          console.log(xex);
      }
  }
  
  for (var i = 0; i < 1000; ++i)
  {
      var xex = new Decimal(-0.3678794411710499+Math.exp(Math.random()*100));
      var x = Decimal.lambertw(xex);
      if (!Decimal.eq_tolerance(xex, x.mul(Decimal.exp(x))))
      {
          console.log(xex);
      }
  }
  
  for (var i = 0; i < 1000; ++i)
  {
      var a = Decimal.randomDecimalForTesting(Math.random() > 0.5 ? 0 : 1);
      var b = Decimal.randomDecimalForTesting(Math.random() > 0.5 ? 0 : 1);
      if (Math.random() > 0.5) { a = a.recip(); }
      if (Math.random() > 0.5) { b = b.recip(); }
      var c = a.add(b).toNumber();
      if (Number.isFinite(c) && !Decimal.eq_tolerance(c, a.toNumber()+b.toNumber()))
      {
          console.log(a + ", " + b);
      }
  }
  
  for (var i = 0; i < 100; ++i)
  {
      var a = Decimal.randomDecimalForTesting(Math.round(Math.random()*4));
      var b = Decimal.randomDecimalForTesting(Math.round(Math.random()*4));
      if (Math.random() > 0.5) { a = a.recip(); }
      if (Math.random() > 0.5) { b = b.recip(); }
      var c = a.mul(b).toNumber();
      if (Number.isFinite(c) && Number.isFinite(a.toNumber()) && Number.isFinite(b.toNumber()) && a.toNumber() != 0 && b.toNumber() != 0 && c != 0 && !Decimal.eq_tolerance(c, a.toNumber()*b.toNumber()))
      {
          console.log("Test 1: " + a + ", " + b);
      }
      else if (!Decimal.mul(a.recip(), b.recip()).eq_tolerance(Decimal.mul(a, b).recip()))
      {
          console.log("Test 3: " + a + ", " + b);
      }
  }
  
  for (var i = 0; i < 10; ++i)
  {
      var a = Decimal.randomDecimalForTesting(Math.round(Math.random()*4));
      var b = Decimal.randomDecimalForTesting(Math.round(Math.random()*4));
      if (Math.random() > 0.5 && a.sign !== 0) { a = a.recip(); }
      if (Math.random() > 0.5 && b.sign !== 0) { b = b.recip(); }
      var c = a.pow(b);
      var d = a.root(b.recip());
      var e = a.pow(b.recip());
      var f = a.root(b);
      
      if (!c.eq_tolerance(d) && a.sign !== 0 && b.sign !== 0)
      {
        console.log("Test 1: " + a + ", " + b);
      }
      if (!e.eq_tolerance(f) && a.sign !== 0 && b.sign !== 0)
      {
        console.log("Test 2: " + a + ", " + b);
      }
  }
  
  for (var i = 0; i < 10; ++i)
  {
      var a = Math.round(Math.random()*18-9);
      var b = Math.round(Math.random()*100-50);
      var c = Math.round(Math.random()*18-9);
      var d = Math.round(Math.random()*100-50);
      console.log("Decimal.pow(Decimal.fromMantissaExponent(" + a + ", " + b + "), Decimal.fromMantissaExponent(" + c + ", " + d + ")).toString()");
  }
  
  */
      
      //Pentation/pentate: The result of tetrating 'height' times in a row. An absurdly strong operator - Decimal.pentate(2, 4.28) and Decimal.pentate(10, 2.37) are already too huge for break_eternity.js!
      // https://en.wikipedia.org/wiki/Pentation
      Decimal.prototype.pentate = function(height = 2, payload = FC_NN(1, 0, 1)) {
        payload = D(payload);
        var oldheight = height;
        height = Math.trunc(height);
        var fracheight = oldheight-height;
        
        //I have no idea if this is a meaningful approximation for pentation to continuous heights, but it is monotonic and continuous.
        if (fracheight !== 0)
        {
          if (payload.eq(Decimal.dOne))
          {
            ++height;
            payload = new Decimal(fracheight);
          }
          else
          {
            if (this.eq(10))
            {
              payload = payload.layeradd10(fracheight);
            }
            else
            {
              payload = payload.layeradd(fracheight, this);
            }
          }
        }
        
        for (var i = 0; i < height; ++i)
        {
          payload = this.tetrate(payload);
          //bail if we're NaN
          if (!isFinite(payload.layer) || !isFinite(payload.mag)) { return payload; }
          //give up after 10 iterations if nothing is happening
          if (i > 10) { return payload; }
        }
        
        return payload;
      }
      
      // trig functions!
      Decimal.prototype.sin = function () {
        if (this.mag < 0) { return this; }
        if (this.layer === 0) { return D(Math.sin(this.sign*this.mag)); }
        return FC_NN(0, 0, 0);
      };
  
      Decimal.prototype.cos = function () {
        if (this.mag < 0) { return Decimal.dOne; }
        if (this.layer === 0) { return D(Math.cos(this.sign*this.mag)); }
        return FC_NN(0, 0, 0);
      };
  
      Decimal.prototype.tan = function () {
        if (this.mag < 0) { return this; }
        if (this.layer === 0) { return D(Math.tan(this.sign*this.mag)); }
        return FC_NN(0, 0, 0);
      };
  
      Decimal.prototype.asin = function () {
        if (this.mag < 0) { return this; }
        if (this.layer === 0) { return D(Math.asin(this.sign*this.mag)); }
        return FC_NN(Number.NaN, Number.NaN, Number.NaN);
      };
  
      Decimal.prototype.acos = function () {
        if (this.mag < 0) { return D(Math.acos(this.toNumber())); }
        if (this.layer === 0) { return D(Math.acos(this.sign*this.mag)); }
        return FC_NN(Number.NaN, Number.NaN, Number.NaN);
      };
  
      Decimal.prototype.atan = function () {
        if (this.mag < 0) { return this; }
        if (this.layer === 0) { return D(Math.atan(this.sign*this.mag)); }
        return D(Math.atan(this.sign*1.8e308));
      };
  
      Decimal.prototype.sinh = function () {
        return this.exp().sub(this.negate().exp()).div(2);
      };
  
      Decimal.prototype.cosh = function () {
        return this.exp().add(this.negate().exp()).div(2);
      };
  
      Decimal.prototype.tanh = function () {
        return this.sinh().div(this.cosh());
      };
  
      Decimal.prototype.asinh = function () {
        return Decimal.ln(this.add(this.sqr().add(1).sqrt()));
      };
  
      Decimal.prototype.acosh = function () {
        return Decimal.ln(this.add(this.sqr().sub(1).sqrt()));
      };
  
      Decimal.prototype.atanh = function () {
        if (this.abs().gte(1)) {
          return FC_NN(Number.NaN, Number.NaN, Number.NaN);
        }
  
        return Decimal.ln(this.add(1).div(D(1).sub(this))).div(2);
      };
      
      /**
       * Joke function from Realm Grinder
       */
      Decimal.prototype.ascensionPenalty = function (ascensions) {
        if (ascensions === 0) {
          return this;
        }
  
        return this.root(Decimal.pow(10, ascensions));
      };
      
      /**
       * Joke function from Cookie Clicker. It's 'egg'
       */
      Decimal.prototype.egg = function () {
        return this.add(9);
      };
      
      Decimal.prototype.lessThanOrEqualTo = function (other) {
        return this.cmp(other) < 1;
      };
  
      Decimal.prototype.lessThan = function (other) {
        return this.cmp(other) < 0;
      };
  
      Decimal.prototype.greaterThanOrEqualTo = function (other) {
        return this.cmp(other) > -1;
      };
  
      Decimal.prototype.greaterThan = function (other) {
        return this.cmp(other) > 0;
      };
  
      return Decimal;
    }();
  
      Decimal.dZero = FC_NN(0, 0, 0);
      Decimal.dOne = FC_NN(1, 0, 1);
      Decimal.dNegOne = FC_NN(-1, 0, 1);
      Decimal.dTwo = FC_NN(1, 0, 2);
      Decimal.dTen = FC_NN(1, 0, 10);
      Decimal.dNaN = FC_NN(Number.NaN, Number.NaN, Number.NaN);
      Decimal.dInf = FC_NN(1, Number.POSITIVE_INFINITY, Number.POSITIVE_INFINITY);
      Decimal.dNegInf = FC_NN(-1, Number.NEGATIVE_INFINITY, Number.NEGATIVE_INFINITY);
    Decimal.dNumberMax = FC(1, 0, Number.MAX_VALUE);
    Decimal.dNumberMin = FC(1, 0, Number.MIN_VALUE);
    
    return Decimal;
  
  }));
  
function exponentialFormat(num, precision, mantissa = true) {commaFormat
    let e = num.log10().floor()
    let m = num.div(Decimal.pow(10, e))
    if (m.toStringWithDecimalPlaces(precision) == 10) {
        m = new Decimal(1)
        e = e.add(1)
    }
    e = (e.gte(1e9) ? format(e, 3) : (e.gte(10000) ? (e, 0) : e.toStringWithDecimalPlaces(0)))
    if (mantissa)
        return m.toStringWithDecimalPlaces(precision) + "e" + e
    else return "e" + e
}

function commaFormat(num, precision) {
    if (num === null || num === undefined) return "NaN"
    if (num.mag < 0.001) return (0).toFixed(precision)
    let init = num.toStringWithDecimalPlaces(precision)
    let portions = init.split(".")
    portions[0] = portions[0].replace(/(\d)(?=(\d\d\d)+(?!\d))/g, "$1,")
    if (portions.length == 1) return portions[0]
    return portions[0] + "." + portions[1]
}


function regularFormat(num, precision) {
    if (num === null || num === undefined) return "NaN"
    if (num.mag < 0.0001) return (0).toFixed(precision)
    if (num.mag < 0.1 && precision !==0) precision = Math.max(precision, 4)
    return num.toStringWithDecimalPlaces(precision)
}

function fixValue(x, y = 0) {
    return x || new Decimal(y)
}

function sumValues(x) {
    x = Object.values(x)
    if (!x[0]) return decimalZero
    return x.reduce((a, b) => Decimal.add(a, b))
}

function format(decimal, precision = 2, small) {
    small = small
    decimal = new Decimal(decimal)
    if (decimal.sign < 0) return "-" + format(decimal.neg(), precision, small)
    if (decimal.mag == Number.POSITIVE_INFINITY) return "Infinity"
    if (decimal.gte("eeee1000")) {
        var slog = decimal.slog()
        return Decimal.pow(10, slog.sub(slog.floor())).toStringWithDecimalPlaces(3) + "F" + commaFormat(slog.floor(), 0)
    }
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(54)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(54)).mul(7)))),precision)+"那摩怛罗"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(53)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(53)).mul(7)))),precision)+"伽摩怛罗"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(52)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(52)).mul(7)))),precision)+"勃摩怛罗"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(51)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(51)).mul(7)))),precision)+"阿摩怛罗"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(50)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(50)).mul(7)))),precision)+"极量"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(49)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(49)).mul(7)))),precision)+"不动"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(48)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(48)).mul(7)))),precision)+"离娇慢"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(47)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(47)).mul(7)))),precision)+"摩鲁摩"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(46)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(46)).mul(7)))),precision)+"翳罗陀"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(45)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(45)).mul(7)))),precision)+"忏慕陀"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(44)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(44)).mul(7)))),precision)+"摩鲁陀"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(43)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(43)).mul(7)))),precision)+"诃鲁那"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(42)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(42)).mul(7)))),precision)+"达罗步陀"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(41)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(41)).mul(7)))),precision)+"奚鲁伽"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(40)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(40)).mul(7)))),precision)+"诃理三"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(39)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(39)).mul(7)))),precision)+"诃理蒲"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(38)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(38)).mul(7)))),precision)+"一动"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(37)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(37)).mul(7)))),precision)+"诃理婆"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(36)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(36)).mul(7)))),precision)+"泥罗婆"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(35)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(35)).mul(7)))),precision)+"最妙"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(34)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(34)).mul(7)))),precision)+"高出"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(33)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(33)).mul(7)))),precision)+"周广"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(32)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(32)).mul(7)))),precision)+"伺察"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(31)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(31)).mul(7)))),precision)+"奚婆罗"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(30)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(30)).mul(7)))),precision)+"毗睹罗"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(29)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(29)).mul(7)))),precision)+"三末耶"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(28)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(28)).mul(7)))),precision)+"颠倒"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(27)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(27)).mul(7)))),precision)+"异路"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(26)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(26)).mul(7)))),precision)+"一持"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(25)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(25)).mul(7)))),precision)+"称量"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(24)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(24)).mul(7)))),precision)+"毗佉担"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(23)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(23)).mul(7)))),precision)+"毗薄底"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(22)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(22)).mul(7)))),precision)+"毗婆诃"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(21)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(21)).mul(7)))),precision)+"毗素陀"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(20)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(20)).mul(7)))),precision)+"毗盛伽"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(19)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(19)).mul(7)))),precision)+"毗赡婆"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(18)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(18)).mul(7)))),precision)+"毗萨罗"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(17)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(17)).mul(7)))),precision)+"僧羯罗摩"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(16)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(16)).mul(7)))),precision)+"毗伽婆"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(15)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(15)).mul(7)))),precision)+"毗罗伽"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(14)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(14)).mul(7)))),precision)+"弥伽婆"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(13)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(13)).mul(7)))),precision)+"阿婆钤"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(12)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(12)).mul(7)))),precision)+"祢摩"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(11)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(11)).mul(7)))),precision)+"普摩"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(10)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(10)).mul(7)))),precision)+"界分"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(9)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(9)).mul(7)))),precision)+"多婆罗"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(8)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(8)).mul(7)))),precision)+"阿婆罗"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(7)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(7)).mul(7)))),precision)+"摩婆罗"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(6)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(6)).mul(7)))),precision)+"最胜"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(5)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(5)).mul(7)))),precision)+"阿伽罗"
    else if (decimal.gte(new Decimal(10).pow((new Decimal(2).pow(4)).mul(7)))) return exponentialFormat((decimal.div(new Decimal(10).pow((new Decimal(2).pow(4)).mul(7)))),precision)+"矜羯罗"
    else if (decimal.gte(1e72)) return exponentialFormat((decimal.div(1e68)),precision)+"无量大数"
    else if (decimal.gte(1e68)) return formatWhole((Math.floor(decimal / 1e68)))+"无量大数"+formatWhole((Math.floor(decimal / 1e64 % 10000)))+"不可思议"+formatWhole((Math.floor(decimal / 1e60 % 10000)))+"那由他"
    else if (decimal.gte(1e64)) return formatWhole((Math.floor(decimal / 1e64)))+"不可思议"+formatWhole((Math.floor(decimal / 1e60 % 10000)))+"那由他"+formatWhole((Math.floor(decimal / 1e56 % 10000)))+"阿僧祇"
    else if (decimal.gte(1e60)) return formatWhole((Math.floor(decimal / 1e60)))+"那由他"+formatWhole((Math.floor(decimal / 1e56 % 10000)))+"阿僧祇"+formatWhole((Math.floor(decimal / 1e52 % 10000)))+"恒河沙"
    else if (decimal.gte(1e56)) return formatWhole((Math.floor(decimal / 1e56)))+"阿僧祇"+formatWhole((Math.floor(decimal / 1e52 % 10000)))+"恒河沙"+formatWhole((Math.floor(decimal / 1e48 % 10000)))+"极"
    else if (decimal.gte(1e52)) return formatWhole((Math.floor(decimal / 1e52)))+"恒河沙"+formatWhole((Math.floor(decimal / 1e48 % 10000)))+"极"+formatWhole((Math.floor(decimal / 1e44 % 10000)))+"载"
    else if (decimal.gte(1e48)) return formatWhole((Math.floor(decimal / 1e48)))+"极"+formatWhole((Math.floor(decimal / 1e44 % 10000)))+"载"+formatWhole((Math.floor(decimal / 1e40 % 10000)))+"正"
    else if (decimal.gte(1e44)) return formatWhole((Math.floor(decimal / 1e44)))+"载"+formatWhole((Math.floor(decimal / 1e40 % 10000)))+"正"+formatWhole((Math.floor(decimal / 1e36 % 10000)))+"涧"
    else if (decimal.gte(1e40)) return formatWhole((Math.floor(decimal / 1e40)))+"正"+formatWhole((Math.floor(decimal / 1e36 % 10000)))+"涧"+formatWhole((Math.floor(decimal / 1e32 % 10000)))+"沟"
    else if (decimal.gte(1e36)) return formatWhole((Math.floor(decimal / 1e36)))+"涧"+formatWhole((Math.floor(decimal / 1e32 % 10000)))+"沟"+formatWhole((Math.floor(decimal / 1e28 % 10000)))+"穰"
    else if (decimal.gte(1e32)) return formatWhole((Math.floor(decimal / 1e32)))+"沟"+formatWhole((Math.floor(decimal / 1e28 % 10000)))+"穰"+formatWhole((Math.floor(decimal / 1e24 % 10000)))+"杼"
    else if (decimal.gte(1e28)) return formatWhole((Math.floor(decimal / 1e28)))+"穰"+formatWhole((Math.floor(decimal / 1e24 % 10000)))+"杼"+formatWhole((Math.floor(decimal / 1e20 % 10000)))+"垓"
    else if (decimal.gte(1e24)) return formatWhole((Math.floor(decimal / 1e24)))+"杼"+formatWhole((Math.floor(decimal / 1e20 % 10000)))+"垓"+formatWhole((Math.floor(decimal / 1e16 % 10000)))+"京"
    else if (decimal.gte(1e20)) return formatWhole((Math.floor(decimal / 1e20)))+"垓"+formatWhole((Math.floor(decimal / 1e16 % 10000)))+"京"+formatWhole((Math.floor(decimal / 1e12 % 10000)))+"兆"
    else if (decimal.gte(1e16)) return formatWhole((Math.floor(decimal / 1e16)))+"京"+formatWhole((Math.floor(decimal / 1e12 % 10000)))+"兆"+formatWhole((Math.floor(decimal / 1e8 % 10000)))+"亿"
    else if (decimal.gte(1e12)) return formatWhole((Math.floor(decimal / 1e12)))+"兆"+formatWhole((Math.floor(decimal / 1e8 % 10000)))+"亿"+formatWhole((Math.floor(decimal / 1e4 % 10000)))+"万"
    else if (decimal.gte(1e8)) return formatWhole((Math.floor(decimal / 1e8)))+"亿"+formatWhole((Math.floor(decimal / 1e4 % 10000)))+"万"+formatWhole((Math.floor(decimal / 1 % 10000)))+""
    else if (decimal.gte(1e4)) return formatWhole((Math.floor(decimal / 1e4)))+"万"+formatWhole((Math.floor(decimal / 1 % 10000)))+""
    else if (decimal.gte(0.0001) || !small) return regularFormat(decimal, precision)
    else if (decimal.eq(0)) return (0).toFixed(precision)

    decimal = invertOOM(decimal)
    let val = ""
    if (decimal.lt("1e1000")){
        val = exponentialFormat(decimal, precision)
        return val.replace(/([^(?:e|F)]*)$/, '-$1')
    }
    else   
        return format(decimal, precision) + "⁻¹"

}

function formatWhole(decimal) {
    decimal = new Decimal(decimal)
    if (decimal.gte(1e9)) return format(decimal, 2)
    if (decimal.lte(0.99) && !decimal.eq(0)) return format(decimal, 2)
    return format(decimal, 0)
}

function formatTime(s) {
    if (s < 60) return format(s) + "s"
    else if (s < 3600) return formatWhole(Math.floor(s / 60)) + "m " + format(s % 60) + "s"
    else if (s < 86400) return formatWhole(Math.floor(s / 3600)) + "h " + formatWhole(Math.floor(s / 60) % 60) + "m " + format(s % 60) + "s"
    else if (s < 31536000) return formatWhole(Math.floor(s / 86400) % 365) + "d " + formatWhole(Math.floor(s / 3600) % 24) + "h " + formatWhole(Math.floor(s / 60) % 60) + "m " + format(s % 60) + "s"
    else return formatWhole(Math.floor(s / 31536000)) + "y " + formatWhole(Math.floor(s / 86400) % 365) + "d " + formatWhole(Math.floor(s / 3600) % 24) + "h " + formatWhole(Math.floor(s / 60) % 60) + "m " + format(s % 60) + "s"
}

function toPlaces(x, precision, maxAccepted) {
    x = new Decimal(x)
    let result = x.toStringWithDecimalPlaces(precision)
    if (new Decimal(result).gte(maxAccepted)) {
        result = new Decimal(maxAccepted - Math.pow(0.1, precision)).toStringWithDecimalPlaces(precision)
    }
    return result
}

// Will also display very small numbers
function formatSmall(x, precision=2) { 
    return format(x, precision, true)    
}

function invertOOM(x){
    let e = x.log10().ceil()
    let m = x.div(Decimal.pow(10, e))
    e = e.neg()
    x = new Decimal(10).pow(e).times(m)

    return x
}
window.onload=function(){
  option1()
if(userData.Lingqi === undefined) userData.Lingqi = new Decimal(0)
if(userData.totalLingqi === undefined) userData.totalLingqi = new Decimal(0)
if(userData.LingqiPerSec === undefined) userData.LingqiPerSec = new Decimal(0)
if(userData.XiuweiPerSec === undefined) userData.XiuweiPerSec = new Decimal(0)
if(userData.Lingshi === undefined) userData.Lingshi = new Decimal(0)
if(userData.Xiuwei === undefined) userData.Xiuwei = new Decimal(0)
if(userData.Jingjieid === undefined) userData.Jingjieid = new Decimal(1)
if(userData.Gongfa1Level === undefined) userData.Gongfa1Level = new Decimal(0)
if(userData.Gongfa2Level === undefined) userData.Gongfa2Level = new Decimal(0)
if(userData.Gongfa3Level === undefined) userData.Gongfa3Level = new Decimal(0)
if(userData.Gongfa4Level === undefined) userData.Gongfa4Level = new Decimal(0)
if(userData.Gongfa5Level === undefined) userData.Gongfa5Level = new Decimal(0)
if(userData.Gongfa6Level === undefined) userData.Gongfa6Level = new Decimal(0)
if(userData.Gongfa7Level === undefined) userData.Gongfa7Level = new Decimal(0)
if(userData.Gongfa8Level === undefined) userData.Gongfa8Level = new Decimal(0)
if(userData.Gongfa9Level === undefined) userData.Gongfa9Level = new Decimal(0)
if(userData.Julingzhen1Level === undefined) userData.Julingzhen1Level = new Decimal(0)
if(userData.Julingzhen1Freeze === undefined) userData.Julingzhen1Freeze = new Decimal(3)
if(userData.Julingzhen1FreezeLimit === undefined) userData.Julingzhen1FreezeLimit = new Decimal(3)
if(userData.Julingzhen2Level === undefined) userData.Julingzhen2Level = new Decimal(0)
if(userData.Julingzhen2Freeze === undefined) userData.Julingzhen2Freeze = new Decimal(6)
if(userData.Julingzhen2FreezeLimit === undefined) userData.Julingzhen2FreezeLimit = new Decimal(6)
if(userData.Julingzhen3Level === undefined) userData.Julingzhen3Level = new Decimal(0)
if(userData.Julingzhen3Freeze === undefined) userData.Julingzhen3Freeze = new Decimal(10)
if(userData.Julingzhen3FreezeLimit === undefined) userData.Julingzhen3FreezeLimit = new Decimal(10)
if(userData.Julingzhen4Level === undefined) userData.Julingzhen4Level = new Decimal(0)
if(userData.Julingzhen4Freeze === undefined) userData.Julingzhen4Freeze = new Decimal(15)
if(userData.Julingzhen4FreezeLimit === undefined) userData.Julingzhen4FreezeLimit = new Decimal(15)
if(userData.Julingzhen5Level === undefined) userData.Julingzhen5Level = new Decimal(0)
if(userData.Julingzhen5Freeze === undefined) userData.Julingzhen5Freeze = new Decimal(35)
if(userData.Julingzhen5FreezeLimit === undefined) userData.Julingzhen5FreezeLimit = new Decimal(35)
if(userData.Julingzhen6Level === undefined) userData.Julingzhen6Level = new Decimal(0)
if(userData.Julingzhen6Freeze === undefined) userData.Julingzhen6Freeze = new Decimal(50)
if(userData.Julingzhen6FreezeLimit === undefined) userData.Julingzhen6FreezeLimit = new Decimal(50)
if(userData.Julingzhen7Level === undefined) userData.Julingzhen7Level = new Decimal(0)
if(userData.Julingzhen7Freeze === undefined) userData.Julingzhen7Freeze = new Decimal(75)
if(userData.Julingzhen7FreezeLimit === undefined) userData.Julingzhen7FreezeLimit = new Decimal(75)
if(userData.Julingzhen8Level === undefined) userData.Julingzhen8Level = new Decimal(0)
if(userData.Julingzhen8Freeze === undefined) userData.Julingzhen8Freeze = new Decimal(100)
if(userData.Julingzhen8FreezeLimit === undefined) userData.Julingzhen8FreezeLimit = new Decimal(100)
if(userData.Jingjie === undefined) userData.Jingjie = '练气1阶'
if(userData.XiulianBalance === undefined) userData.XiulianBalance = new Decimal(0)
getAchievement()
}
var Lingqi;
var LingqiPerSec;
var Xiuwei;
var XiuweiPerSec;
var Lingshi;
var Jingjieid;

Gongfa1cost = [new Decimal(100),new Decimal(500),new Decimal(1000),new Decimal(3000),new Decimal(8000),new Decimal(15000),new Decimal(30000),new Decimal(60000),new Decimal(100000),new Decimal(150000),new Decimal("eee100")]
Gongfa2cost = [new Decimal(5000),new Decimal(9000),new Decimal(20000),new Decimal(35000),new Decimal(100000),new Decimal(250000),new Decimal(400000),new Decimal(1000000),new Decimal(1550000),new Decimal(2500000),new Decimal("eee100")]
Gongfa3cost = [new Decimal(30000),new Decimal(50000),new Decimal(80000),new Decimal(150000),new Decimal(300000),new Decimal(500000),new Decimal(800000),new Decimal(2000000),new Decimal(3100000),new Decimal(5000000),new Decimal("eee100")]
Gongfa4cost = [new Decimal(80000),new Decimal(120000),new Decimal(200000),new Decimal(300000),new Decimal(650000),new Decimal(1250000),new Decimal(2700000),new Decimal(5000000),new Decimal(9000000),new Decimal(15000000),new Decimal("eee100")]
Gongfa5cost = [new Decimal(1000000),new Decimal(1500000),new Decimal(2100000),new Decimal(3000000),new Decimal(3800000),new Decimal(5500000),new Decimal(8000000),new Decimal(13000000),new Decimal(20000000),new Decimal(35000000),new Decimal("eee100")]
Gongfa6cost = [new Decimal(5000000),new Decimal(6500000),new Decimal(9000000),new Decimal(12000000),new Decimal(18000000),new Decimal(35000000),new Decimal(70000000),new Decimal(130000000),new Decimal(190000000),new Decimal(250000000),new Decimal("eee100")]
Gongfa7cost = [new Decimal(30000000),new Decimal(50000000),new Decimal(75000000),new Decimal(100000000),new Decimal(135000000),new Decimal(180000000),new Decimal(250000000),new Decimal(400000000),new Decimal(800000000),new Decimal(1000000000),new Decimal("eee100")]
Gongfa8cost = [new Decimal(8e7),new Decimal(1e8),new Decimal(1.5e8),new Decimal(2.2e8),new Decimal(3.3e8),new Decimal(4.5e8),new Decimal(6e8),new Decimal(8e8),new Decimal(1.2e9),new Decimal(2e9),new Decimal("eee100")]
Gongfa9cost = [new Decimal(1.5e8),new Decimal(2e8),new Decimal(3e8),new Decimal(4.5e8),new Decimal(6.8e8),new Decimal(9e8),new Decimal(1.2e9),new Decimal(1.6e9),new Decimal(3e9),new Decimal(5e9),new Decimal("eee100")]
Julingzhen1cost = [new Decimal(1000),new Decimal(3000),new Decimal(8000),new Decimal(15000),new Decimal(30000),new Decimal(60000),new Decimal(100000),new Decimal(150000),new Decimal(220000),new Decimal(300000),new Decimal(400000),new Decimal(550000),new Decimal(750000),new Decimal(1100000),new Decimal(1500000),new Decimal(2000000),new Decimal(3000000),new Decimal(4500000),new Decimal(7000000),new Decimal("eee100")]
Julingzhen2cost = [new Decimal(5000),new Decimal(10000),new Decimal(30000),new Decimal(100000),new Decimal(250000),new Decimal(400000),new Decimal(600000),new Decimal(800000),new Decimal(1000000),new Decimal(1300000),new Decimal(1750000),new Decimal(2300000),new Decimal(3000000),new Decimal(4000000),new Decimal(5500000),new Decimal(7500000),new Decimal(1e7),new Decimal(1.5e7),new Decimal(2.1e7),new Decimal("eee100")]
Julingzhen3cost = [new Decimal(30000),new Decimal(100000),new Decimal(350000),new Decimal(600000),new Decimal(900000),new Decimal(1300000),new Decimal(1900000),new Decimal(2400000),new Decimal(3000000),new Decimal(3800000),new Decimal(4700000),new Decimal(6000000),new Decimal(8000000),new Decimal(11000000),new Decimal(15000000),new Decimal(20000000),new Decimal(2.7e7),new Decimal(3.5e7),new Decimal(5e7),new Decimal("eee100")]
Julingzhen4cost = [new Decimal(200000),new Decimal(500000),new Decimal(1000000),new Decimal(1700000),new Decimal(2800000),new Decimal(4300000),new Decimal(6000000),new Decimal(8400000),new Decimal(1e7),new Decimal(1.2e7),new Decimal(1.5e7),new Decimal(1.9e7),new Decimal(2.5e7),new Decimal(3.3e7),new Decimal(4e7),new Decimal(5e7),new Decimal(6.5e7),new Decimal(8.5e7),new Decimal(1.1e8),new Decimal("eee100")]
Julingzhen5cost = [new Decimal(1000000),new Decimal(3000000),new Decimal(8000000),new Decimal(1.5e7),new Decimal(3e7),new Decimal(6e7),new Decimal(1e8),new Decimal(1.5e8),new Decimal(2.2e8),new Decimal(3e8),new Decimal(4e8),new Decimal(5.5e8),new Decimal(7.5e8),new Decimal(1e9),new Decimal(1.3e9),new Decimal(1.7e9),new Decimal(2.2e9),new Decimal(2.8e9),new Decimal(3.5e9),new Decimal("eee100")]
Julingzhen6cost = [new Decimal(5000000),new Decimal(1e7),new Decimal(3e7),new Decimal(1e8),new Decimal(2.5e8),new Decimal(4e8),new Decimal(6e8),new Decimal(8e8),new Decimal(1e9),new Decimal(1.3e9),new Decimal(1.75e9),new Decimal(2.3e9),new Decimal(3e9),new Decimal(3.8e9),new Decimal(5e9),new Decimal(6.5e9),new Decimal(8.5e9),new Decimal(1.1e10),new Decimal(1.5e10),new Decimal("eee100")]
Julingzhen7cost = [new Decimal(3e7),new Decimal(1e8),new Decimal(3.5e8),new Decimal(6e8),new Decimal(9e8),new Decimal(1.3e9),new Decimal(1.9e9),new Decimal(2.4e9),new Decimal(3e9),new Decimal(3.8e9),new Decimal(4.6e9),new Decimal(5.5e9),new Decimal(7e9),new Decimal(9e9),new Decimal(1.2e10),new Decimal(1.7e10),new Decimal(2.5e10),new Decimal(3.5e10),new Decimal(5.5e10),new Decimal("eee100")]
Julingzhen8cost = [new Decimal(1e8),new Decimal(5e8),new Decimal(1e9),new Decimal(1.7e9),new Decimal(2.8e9),new Decimal(4.3e9),new Decimal(6e9),new Decimal(8.4e9),new Decimal(1e10),new Decimal(1.2e10),new Decimal(1.5e10),new Decimal(2e10),new Decimal(2.75e10),new Decimal(3.7e10),new Decimal(5e10),new Decimal(7e10),new Decimal(1e11),new Decimal(1.8e11),new Decimal(3e11),new Decimal("eee100")]
Gongfa1effect = [new Decimal(0),new Decimal(20),new Decimal(30),new Decimal(50),new Decimal(100),new Decimal(170),new Decimal(260),new Decimal(420),new Decimal(640),new Decimal(900),new Decimal(1300)]
Gongfa2effect = [new Decimal(0),new Decimal(60),new Decimal(100),new Decimal(160),new Decimal(280),new Decimal(440),new Decimal(700),new Decimal(1100),new Decimal(1600),new Decimal(2400),new Decimal(3600)]
Gongfa3effect = [new Decimal(0),new Decimal(80),new Decimal(180),new Decimal(260),new Decimal(600),new Decimal(1000),new Decimal(1500),new Decimal(2200),new Decimal(3100),new Decimal(4500),new Decimal(6500)]
Gongfa4effect = [new Decimal(0),new Decimal(100),new Decimal(240),new Decimal(300),new Decimal(800),new Decimal(1300),new Decimal(1900),new Decimal(3100),new Decimal(5500),new Decimal(8300),new Decimal(12000)]
Gongfa5effect = [new Decimal(0),new Decimal(1000),new Decimal(1500),new Decimal(2500),new Decimal(5000),new Decimal(8500),new Decimal(13000),new Decimal(21000),new Decimal(32000),new Decimal(45000),new Decimal(65000)]
Gongfa6effect = [new Decimal(0),new Decimal(3000),new Decimal(5000),new Decimal(8000),new Decimal(14000),new Decimal(22000),new Decimal(35000),new Decimal(55000),new Decimal(80000),new Decimal(120000),new Decimal(180000)]
Gongfa7effect = [new Decimal(0),new Decimal(5000),new Decimal(10000),new Decimal(18000),new Decimal(40000),new Decimal(75000),new Decimal(130000),new Decimal(200000),new Decimal(310000),new Decimal(450000),new Decimal(650000)]
Gongfa8effect = [new Decimal(0),new Decimal(8500),new Decimal(16000),new Decimal(30000),new Decimal(60000),new Decimal(110000),new Decimal(190000),new Decimal(310000),new Decimal(550000),new Decimal(830000),new Decimal(1200000)]
Gongfa9effect = [new Decimal(0),new Decimal(15000),new Decimal(25000),new Decimal(45000),new Decimal(70000),new Decimal(150000),new Decimal(250000),new Decimal(380000),new Decimal(650000),new Decimal(990000),new Decimal(1919810)]
Julingzhen1effect = [new Decimal(0),new Decimal(10),new Decimal(20),new Decimal(50),new Decimal(100),new Decimal(170),new Decimal(270),new Decimal(440),new Decimal(700),new Decimal(1000),new Decimal(1400),new Decimal(1800),new Decimal(2300),new Decimal(2900),new Decimal(3700),new Decimal(4500),new Decimal(5500),new Decimal(6600),new Decimal(8000),new Decimal(9500)]
Julingzhen2effect = [new Decimal(0),new Decimal(50),new Decimal(80),new Decimal(130),new Decimal(210),new Decimal(360),new Decimal(630),new Decimal(980),new Decimal(1470),new Decimal(2100),new Decimal(2900),new Decimal(3900),new Decimal(5150),new Decimal(6700),new Decimal(8700),new Decimal(11300),new Decimal(14400),new Decimal(18200),new Decimal(22700),new Decimal(30000)]
Julingzhen3effect = [new Decimal(0),new Decimal(100),new Decimal(150),new Decimal(250),new Decimal(420),new Decimal(680),new Decimal(1050),new Decimal(1540),new Decimal(2170),new Decimal(3000),new Decimal(4050),new Decimal(5550),new Decimal(7550),new Decimal(10250),new Decimal(13800),new Decimal(18600),new Decimal(25100),new Decimal(33600),new Decimal(44600),new Decimal(59000)]
Julingzhen4effect = [new Decimal(0),new Decimal(250),new Decimal(350),new Decimal(600),new Decimal(900),new Decimal(1350),new Decimal(2000),new Decimal(2900),new Decimal(4500),new Decimal(6500),new Decimal(9000),new Decimal(12300),new Decimal(16300),new Decimal(21300),new Decimal(27600),new Decimal(35600),new Decimal(46000),new Decimal(59000),new Decimal(75500),new Decimal(96000)]
Julingzhen5effect = [new Decimal(0),new Decimal(1000),new Decimal(2000),new Decimal(5000),new Decimal(10000),new Decimal(17000),new Decimal(27000),new Decimal(44000),new Decimal(70000),new Decimal(90000),new Decimal(140000),new Decimal(213000),new Decimal(313000),new Decimal(450000),new Decimal(630000),new Decimal(850000),new Decimal(1130000),new Decimal(1480000),new Decimal(1910000),new Decimal(2300000)]
Julingzhen6effect = [new Decimal(0),new Decimal(5000),new Decimal(8000),new Decimal(13000),new Decimal(21000),new Decimal(36000),new Decimal(63000),new Decimal(98000),new Decimal(147000),new Decimal(210000),new Decimal(285000),new Decimal(375000),new Decimal(495000),new Decimal(650000),new Decimal(850000),new Decimal(1000000),new Decimal(1300000),new Decimal(1750000),new Decimal(2300000),new Decimal(3200000)]
Julingzhen7effect = [new Decimal(0),new Decimal(10000),new Decimal(15000),new Decimal(25000),new Decimal(42000),new Decimal(68000),new Decimal(105000),new Decimal(154000),new Decimal(217000),new Decimal(300000),new Decimal(410000),new Decimal(550000),new Decimal(730000),new Decimal(960000),new Decimal(1260000),new Decimal(1650000),new Decimal(2150000),new Decimal(2800000),new Decimal(3650000),new Decimal(4800000)]
Julingzhen8effect = [new Decimal(0),new Decimal(25000),new Decimal(35000),new Decimal(60000),new Decimal(90000),new Decimal(135000),new Decimal(200000),new Decimal(290000),new Decimal(450000),new Decimal(650000),new Decimal(920000),new Decimal(1270000),new Decimal(1770000),new Decimal(2470000),new Decimal(3500000),new Decimal(5000000),new Decimal(7200000),new Decimal(10200000),new Decimal(14200000),new Decimal(17000000)]
JingjieLimit = [new Decimal(0),new Decimal(1000),new Decimal(2000),new Decimal(3500),new Decimal(5500),new Decimal(8500),new Decimal(13000),new Decimal(20000),new Decimal(30000),new Decimal(45000),new Decimal(65000),new Decimal(90000),new Decimal(120000),new Decimal(190000),new Decimal(250000),new Decimal(300000),new Decimal(380000),new Decimal(450000),new Decimal(550000),new Decimal(700000),new Decimal(870000),new Decimal(1050000),new Decimal(1400000),new Decimal(2000000),new Decimal(2600000),new Decimal(3500000),new Decimal(4500000),new Decimal(6000000),new Decimal(8000000),new Decimal(11000000),new Decimal(13000000),new Decimal(15000000),new Decimal(18000000),new Decimal(35000000),new Decimal(55000000),new Decimal(120000000),new Decimal(160000000),new Decimal(200000000),new Decimal(250000000),new Decimal(300000000),new Decimal(370000000),new Decimal(420000000),new Decimal(550000000),new Decimal(640000000),new Decimal(750000000),new Decimal(1e12)]
//将元素变量附加到DOM元素
Lingqi = document.getElementById("Lingqi");
Xiuwei = document.getElementById("Xiuwei");
Jingjieid = document.getElementById("Jingjieid");
Lingshi = document.getElementById("Lingshi");
LingqiPerSec = document.getElementById("LingqiPerSec");
var userData = {};
//加载用户数据
if(localStorage.getItem("gameStateData") === null){
    //创建新的游戏数据
    userData = {
    Lingqi: new Decimal(100),
    LingqiPerSec: new Decimal(0),
    Lingshi: new Decimal(0),
    Xiuwei: new Decimal(0),
    XiuweiPerSec: new Decimal(0),
    Jingjieid: new Decimal(1),
    Gongfa1Level: new Decimal(0),
    Gongfa2Level: new Decimal(0),
    Gongfa3Level: new Decimal(0),
    Gongfa4Level: new Decimal(0),
    Gongfa5Level: new Decimal(0),
    Gongfa6Level: new Decimal(0),
    Gongfa7Level: new Decimal(0),
    Gongfa8Level: new Decimal(0),
    Gongfa9Level: new Decimal(0),
    Julingzhen1Level: new Decimal(0),
    Julingzhen1Freeze: new Decimal(3),
    Julingzhen1FreezeLimit: new Decimal(3),
    Julingzhen2Level: new Decimal(0),
    Julingzhen2Freeze: new Decimal(6),
    Julingzhen2FreezeLimit: new Decimal(6),
    Julingzhen3Level: new Decimal(0),
    Julingzhen3Freeze: new Decimal(10),
    Julingzhen3FreezeLimit: new Decimal(10),
    Julingzhen4Level: new Decimal(0),
    Julingzhen4Freeze: new Decimal(15),
    Julingzhen4FreezeLimit: new Decimal(15),
    Julingzhen5Level: new Decimal(0),
    Julingzhen5Freeze: new Decimal(35),
    Julingzhen5FreezeLimit: new Decimal(35),
    Julingzhen6Level: new Decimal(0),
    Julingzhen6Freeze: new Decimal(50),
    Julingzhen6FreezeLimit: new Decimal(50),
    Julingzhen7Level: new Decimal(0),
    Julingzhen7Freeze: new Decimal(75),
    Julingzhen7FreezeLimit: new Decimal(75),
    Julingzhen8Level: new Decimal(0),
    Julingzhen8Freeze: new Decimal(100),
    Julingzhen8FreezeLimit: new Decimal(100),
    }
}
else{
    //加载已经保存的游戏数据
    userData = JSON.parse(window.localStorage["gameStateData"]);
}
//硬重置
function reset(){
    //将已经保存的用户数据恢复到初始状态并且刷新页面
    userData = {
    Lingqi: new Decimal(100),
    totalLingqi: new Decimal(100),
    LingqiPerSec: new Decimal(0),
    XiuweiPerSec: new Decimal(0),
    Lingshi: new Decimal(0),
    Xiuwei: new Decimal(0),
    Jingjieid: new Decimal(1),
    Gongfa1Level: new Decimal(0),
    Gongfa2Level: new Decimal(0),
    Gongfa3Level: new Decimal(0),
    Gongfa4Level: new Decimal(0),
    Gongfa5Level: new Decimal(0),
    Gongfa6Level: new Decimal(0),
    Gongfa7Level: new Decimal(0),
    Gongfa8Level: new Decimal(0),
    Gongfa9Level: new Decimal(0),
    Julingzhen1Level: new Decimal(0),
    Julingzhen1Freeze: new Decimal(3),
    Julingzhen1FreezeLimit: new Decimal(3),
    Julingzhen2Level: new Decimal(0),
    Julingzhen2Freeze: new Decimal(6),
    Julingzhen2FreezeLimit: new Decimal(6),
    Julingzhen3Level: new Decimal(0),
    Julingzhen3Freeze: new Decimal(10),
    Julingzhen3FreezeLimit: new Decimal(10),
    Julingzhen4Level: new Decimal(0),
    Julingzhen4Freeze: new Decimal(15),
    Julingzhen4FreezeLimit: new Decimal(15),
    Julingzhen5Level: new Decimal(0),
    Julingzhen5Freeze: new Decimal(35),
    Julingzhen5FreezeLimit: new Decimal(35),
    Julingzhen6Level: new Decimal(0),
    Julingzhen6Freeze: new Decimal(50),
    Julingzhen6FreezeLimit: new Decimal(50),
    Julingzhen7Level: new Decimal(0),
    Julingzhen7Freeze: new Decimal(75),
    Julingzhen7FreezeLimit: new Decimal(75),
    Julingzhen8Level: new Decimal(0),
    Julingzhen8Freeze: new Decimal(100),
    Julingzhen8FreezeLimit: new Decimal(100),
    Jingjie:'练气1阶',
    XiulianBalance: new Decimal(0),
    }
    localStorage.setItem("gameStateData", JSON.stringify(userData));
    location.reload();
}
//保存游戏
function save(){
    userData = {
        Lingqi: new Decimal(userData.Lingqi),
        totalLingqi: new Decimal(userData.totalLingqi),
        Lingshi: new Decimal(userData.Lingshi),
        Xiuwei: new Decimal(userData.Xiuwei),
        Jingjieid: new Decimal(userData.Jingjieid),
        LingqiPerSec: new Decimal(userData.LingqiPerSec),
        XiuweiPerSec: new Decimal(userData.XiuweiPerSec),
        Gongfa1Level: new Decimal(userData.Gongfa1Level),
        Gongfa2Level: new Decimal(userData.Gongfa2Level),
        Gongfa3Level: new Decimal(userData.Gongfa3Level),
        Gongfa4Level: new Decimal(userData.Gongfa4Level),
        Gongfa5Level: new Decimal(userData.Gongfa5Level),
        Gongfa6Level: new Decimal(userData.Gongfa6Level),
        Gongfa7Level: new Decimal(userData.Gongfa7Level),
        Gongfa8Level: new Decimal(userData.Gongfa8Level),
        Gongfa9Level: new Decimal(userData.Gongfa9Level),
        Julingzhen1Level: new Decimal(userData.Julingzhen1Level),
        Julingzhen1Freeze: new Decimal(userData.Julingzhen1Freeze),
        Julingzhen1FreezeLimit: new Decimal(userData.Julingzhen1FreezeLimit),
        Julingzhen2Level: new Decimal(userData.Julingzhen2Level),
        Julingzhen2Freeze: new Decimal(userData.Julingzhen2Freeze),
        Julingzhen2FreezeLimit: new Decimal(userData.Julingzhen2FreezeLimit),
        Julingzhen3Level: new Decimal(userData.Julingzhen3Level),
        Julingzhen3Freeze: new Decimal(userData.Julingzhen3Freeze),
        Julingzhen3FreezeLimit: new Decimal(userData.Julingzhen3FreezeLimit),
        Julingzhen4Level: new Decimal(userData.Julingzhen4Level),
        Julingzhen4Freeze: new Decimal(userData.Julingzhen4Freeze),
        Julingzhen4FreezeLimit: new Decimal(userData.Julingzhen4FreezeLimit),
        Julingzhen5Level: new Decimal(userData.Julingzhen5Level),
        Julingzhen5Freeze: new Decimal(userData.Julingzhen5Freeze),
        Julingzhen5FreezeLimit: new Decimal(userData.Julingzhen5FreezeLimit),
        Julingzhen6Level: new Decimal(userData.Julingzhen6Level),
        Julingzhen6Freeze: new Decimal(userData.Julingzhen6Freeze),
        Julingzhen6FreezeLimit: new Decimal(userData.Julingzhen6FreezeLimit),
        Julingzhen7Level: new Decimal(userData.Julingzhen7Level),
        Julingzhen7Freeze: new Decimal(userData.Julingzhen7Freeze),
        Julingzhen7FreezeLimit: new Decimal(userData.Julingzhen7FreezeLimit),
        Julingzhen8Level: new Decimal(userData.Julingzhen8Level),
        Julingzhen8Freeze: new Decimal(userData.Julingzhen8Freeze),
        Julingzhen8FreezeLimit: new Decimal(userData.Julingzhen8FreezeLimit),

        Jingjie: userData.Jingjie,
        XiulianBalance: new Decimal(userData.XiulianBalance),
    }
    localStorage.setItem("gameStateData", JSON.stringify(userData));
}
//gameloop 
window.setInterval(() => {
    userData.LingqiPerSec = new Decimal(Gongfa1effect[userData.Gongfa1Level]).add(new Decimal(Gongfa2effect[userData.Gongfa2Level])).add(new Decimal(Gongfa3effect[userData.Gongfa3Level])).add(new Decimal(Gongfa4effect[userData.Gongfa4Level])).add(new Decimal(Gongfa5effect[userData.Gongfa5Level])).add(new Decimal(Gongfa6effect[userData.Gongfa6Level])).add(new Decimal(Gongfa7effect[userData.Gongfa7Level])).add(new Decimal(Gongfa8effect[userData.Gongfa8Level])).add(new Decimal(Gongfa9effect[userData.Gongfa9Level]))
    
    userData.Lingqi = userData.Lingqi.add(new Decimal(userData.LingqiPerSec).mul(0.5))
    userData.totalLingqi = userData.totalLingqi.add(new Decimal(userData.LingqiPerSec).mul(0.5))
    document.getElementById("ShowUpgradeLevel1").innerHTML = format(userData.Gongfa1Level)
    document.getElementById("ShowUpgradePrice1").innerHTML = Gongfa1cost[userData.Gongfa1Level].gte("eee10")?"已臻圆满！":"价格"+format(Gongfa1cost[userData.Gongfa1Level])+" 灵气"
    document.getElementById("ShowUpgradeEffect1").innerHTML = format(Gongfa1effect[userData.Gongfa1Level])
    document.getElementById("ShowUpgradeLevel2").innerHTML = format(userData.Gongfa2Level)
    document.getElementById("ShowUpgradePrice2").innerHTML = Gongfa2cost[userData.Gongfa2Level].gte("eee10")?"已臻圆满！":"价格"+format(Gongfa2cost[userData.Gongfa2Level])+" 灵气"
    document.getElementById("ShowUpgradeEffect2").innerHTML = format(Gongfa2effect[userData.Gongfa2Level])
    document.getElementById("ShowUpgradeLevel3").innerHTML = format(userData.Gongfa3Level)
    document.getElementById("ShowUpgradePrice3").innerHTML = Gongfa3cost[userData.Gongfa3Level].gte("eee10")?"已臻圆满！":"价格"+format(Gongfa3cost[userData.Gongfa3Level])+" 灵气"
    document.getElementById("ShowUpgradeEffect3").innerHTML = format(Gongfa3effect[userData.Gongfa3Level])
    document.getElementById("ShowUpgradeLevel4").innerHTML = format(userData.Gongfa4Level)
    document.getElementById("ShowUpgradePrice4").innerHTML = Gongfa4cost[userData.Gongfa4Level].gte("eee10")?"已臻圆满！":"价格"+format(Gongfa4cost[userData.Gongfa4Level])+" 灵气"
    document.getElementById("ShowUpgradeEffect4").innerHTML = format(Gongfa4effect[userData.Gongfa4Level])
    document.getElementById("ShowUpgradeLevel5").innerHTML = format(userData.Gongfa5Level)
    document.getElementById("ShowUpgradePrice5").innerHTML = Gongfa5cost[userData.Gongfa5Level].gte("eee10")?"已臻圆满！":"价格"+format(Gongfa5cost[userData.Gongfa5Level])+" 灵气"
    document.getElementById("ShowUpgradeEffect5").innerHTML = format(Gongfa5effect[userData.Gongfa5Level])
    document.getElementById("ShowUpgradeLevel6").innerHTML = format(userData.Gongfa6Level)
    document.getElementById("ShowUpgradePrice6").innerHTML = Gongfa6cost[userData.Gongfa6Level].gte("eee10")?"已臻圆满！":"价格"+format(Gongfa6cost[userData.Gongfa6Level])+" 灵气"
    document.getElementById("ShowUpgradeEffect6").innerHTML = format(Gongfa6effect[userData.Gongfa6Level])
    document.getElementById("ShowUpgradeLevel7").innerHTML = format(userData.Gongfa7Level)
    document.getElementById("ShowUpgradePrice7").innerHTML = userData.Gongfa7Level>=(10)?"已臻圆满！":"价格"+format(Gongfa7cost[userData.Gongfa7Level])+" 灵气"
    document.getElementById("ShowUpgradeEffect7").innerHTML = format(Gongfa7effect[userData.Gongfa7Level])
    document.getElementById("ShowUpgradeLevel8").innerHTML = format(userData.Gongfa8Level)
    document.getElementById("ShowUpgradePrice8").innerHTML = userData.Gongfa8Level>=(10)?"已臻圆满！":"价格"+format(Gongfa8cost[userData.Gongfa8Level])+" 灵气"
    document.getElementById("ShowUpgradeEffect8").innerHTML = format(Gongfa8effect[userData.Gongfa8Level])
    document.getElementById("ShowUpgradeLevel9").innerHTML = format(userData.Gongfa9Level)
    document.getElementById("ShowUpgradePrice9").innerHTML = userData.Gongfa9Level>=(10)?"已臻圆满！":"价格"+format(Gongfa9cost[userData.Gongfa9Level])+" 灵气"
    document.getElementById("ShowUpgradeEffect9").innerHTML = format(Gongfa9effect[userData.Gongfa9Level])
    document.getElementById("ShowUpgradeLevel21").innerHTML = format(userData.Julingzhen1Level)
    document.getElementById("ShowUpgradePrice21").innerHTML = userData.Julingzhen1Level>=(19)?"(已满级！)":format(Julingzhen1cost[userData.Julingzhen1Level])
    document.getElementById("ShowUpgradeEffect21").innerHTML = format(Julingzhen1effect[userData.Julingzhen1Level])
    document.getElementById("ShowUpgradeLevel22").innerHTML = format(userData.Julingzhen2Level)
    document.getElementById("ShowUpgradePrice22").innerHTML = userData.Julingzhen2Level>=(19)?"(已满级！)":format(Julingzhen2cost[userData.Julingzhen2Level])
    document.getElementById("ShowUpgradeEffect22").innerHTML = format(Julingzhen2effect[userData.Julingzhen2Level])
    document.getElementById("ShowUpgradeLevel23").innerHTML = format(userData.Julingzhen3Level)
    document.getElementById("ShowUpgradePrice23").innerHTML = userData.Julingzhen3Level>=(19)?"(已满级！)":format(Julingzhen3cost[userData.Julingzhen3Level])
    document.getElementById("ShowUpgradeEffect23").innerHTML = format(Julingzhen3effect[userData.Julingzhen3Level])
    document.getElementById("ShowUpgradeLevel24").innerHTML = format(userData.Julingzhen4Level)
    document.getElementById("ShowUpgradePrice24").innerHTML = userData.Julingzhen4Level>=(19)?"(已满级！)":format(Julingzhen4cost[userData.Julingzhen4Level])
    document.getElementById("ShowUpgradeEffect24").innerHTML = format(Julingzhen4effect[userData.Julingzhen4Level])
    document.getElementById("ShowUpgradeLevel25").innerHTML = format(userData.Julingzhen5Level)
    document.getElementById("ShowUpgradePrice25").innerHTML = userData.Julingzhen5Level>=(19)?"(已满级！)":format(Julingzhen5cost[userData.Julingzhen5Level])
    document.getElementById("ShowUpgradeEffect25").innerHTML = format(Julingzhen5effect[userData.Julingzhen5Level])
    document.getElementById("ShowUpgradeLevel26").innerHTML = format(userData.Julingzhen6Level)
    document.getElementById("ShowUpgradePrice26").innerHTML = userData.Julingzhen6Level>=(19)?"(已满级！)":format(Julingzhen6cost[userData.Julingzhen6Level])
    document.getElementById("ShowUpgradeEffect26").innerHTML = format(Julingzhen6effect[userData.Julingzhen6Level])
    document.getElementById("ShowUpgradeLevel27").innerHTML = format(userData.Julingzhen7Level)
    document.getElementById("ShowUpgradePrice27").innerHTML = userData.Julingzhen7Level>=(19)?"(已满级！)":format(Julingzhen7cost[userData.Julingzhen7Level])
    document.getElementById("ShowUpgradeEffect27").innerHTML = format(Julingzhen7effect[userData.Julingzhen7Level])
    document.getElementById("ShowUpgradeLevel28").innerHTML = format(userData.Julingzhen8Level)
    document.getElementById("ShowUpgradePrice28").innerHTML = userData.Julingzhen8Level>=(19)?"(已满级！)":format(Julingzhen8cost[userData.Julingzhen8Level])
    document.getElementById("ShowUpgradeEffect28").innerHTML = format(Julingzhen8effect[userData.Julingzhen8Level])
    document.getElementById("Lingqi").innerHTML = format(userData.Lingqi)
    document.getElementById("Lingshi").innerHTML = format(userData.Lingshi)
    if(userData.Julingzhen1Level.gte(1))userData.Julingzhen1Freeze = userData.Julingzhen1Freeze.sub(0.05)
    if(userData.Julingzhen1Freeze.lt(0.01))userData.Lingshi = userData.Lingshi.add(Julingzhen1effect[userData.Julingzhen1Level]),userData.Julingzhen1Freeze = userData.Julingzhen1FreezeLimit
    if(userData.Julingzhen2Level.gte(1))userData.Julingzhen2Freeze = userData.Julingzhen2Freeze.sub(0.05)
    if(userData.Julingzhen2Freeze.lt(0.01))userData.Lingshi = userData.Lingshi.add(Julingzhen2effect[userData.Julingzhen2Level]),userData.Julingzhen2Freeze = userData.Julingzhen2FreezeLimit
    if(userData.Julingzhen3Level.gte(1))userData.Julingzhen3Freeze = userData.Julingzhen3Freeze.sub(0.05)
    if(userData.Julingzhen3Freeze.lt(0.01))userData.Lingshi = userData.Lingshi.add(Julingzhen3effect[userData.Julingzhen3Level]),userData.Julingzhen3Freeze = userData.Julingzhen3FreezeLimit
    if(userData.Julingzhen4Level.gte(1))userData.Julingzhen4Freeze = userData.Julingzhen4Freeze.sub(0.05)
    if(userData.Julingzhen4Freeze.lt(0.01))userData.Lingshi = userData.Lingshi.add(Julingzhen4effect[userData.Julingzhen4Level]),userData.Julingzhen4Freeze = userData.Julingzhen4FreezeLimit
    if(userData.Julingzhen5Level.gte(1))userData.Julingzhen5Freeze = userData.Julingzhen5Freeze.sub(0.05)
    if(userData.Julingzhen5Freeze.lt(0.01))userData.Lingshi = userData.Lingshi.add(Julingzhen5effect[userData.Julingzhen5Level]),userData.Julingzhen5Freeze = userData.Julingzhen5FreezeLimit
    if(userData.Julingzhen6Level.gte(1))userData.Julingzhen6Freeze = userData.Julingzhen6Freeze.sub(0.05)
    if(userData.Julingzhen6Freeze.lt(0.01))userData.Lingshi = userData.Lingshi.add(Julingzhen6effect[userData.Julingzhen6Level]),userData.Julingzhen6Freeze = userData.Julingzhen6FreezeLimit
    if(userData.Julingzhen7Level.gte(1))userData.Julingzhen7Freeze = userData.Julingzhen7Freeze.sub(0.05)
    if(userData.Julingzhen7Freeze.lt(0.01))userData.Lingshi = userData.Lingshi.add(Julingzhen7effect[userData.Julingzhen7Level]),userData.Julingzhen7Freeze = userData.Julingzhen7FreezeLimit
    if(userData.Julingzhen8Level.gte(1))userData.Julingzhen8Freeze = userData.Julingzhen8Freeze.sub(0.05)
    if(userData.Julingzhen8Freeze.lt(0.01))userData.Lingshi = userData.Lingshi.add(Julingzhen8effect[userData.Julingzhen8Level]),userData.Julingzhen8Freeze = userData.Julingzhen8FreezeLimit
    if(userData.Jingjieid<=15) userData.Jingjie="练气"+userData.Jingjieid+"阶"
    else if(userData.Jingjieid<30) userData.Jingjie= "筑基"+(userData.Jingjieid%15)+"阶"
    else if(userData.Jingjieid==30) userData.Jingjie= "筑基15阶"
    else if(userData.Jingjieid<45) userData.Jingjie= "开光"+(userData.Jingjieid%15)+"阶"
    else if(userData.Jingjieid==45) userData.Jingjie= "筑基15阶"
    if((userData.Xiuwei).gte(JingjieLimit[userData.Jingjieid])) userData.Xiuwei = userData.Xiuwei.sub(JingjieLimit[userData.Jingjieid]),userData.Jingjieid++,getAchievement()
    if(userData.XiulianBalance.gte(0.05)) userData.XiulianBalance = userData.XiulianBalance.sub(0.05),userData.Xiuwei = userData.Xiuwei.add(userData.XiuweiPerSec)
    if(userData.XiulianBalance.lte(0)) userData.XiuweiPerSec = new Decimal(0)
    document.getElementById("UpdateJingjie").innerHTML = userData.Jingjie
    
    option2Appear()
    option3Appear()
    option4Appear()
    ZhujiAppear()
    KaiguangAppear()
    option3_2Appear()
    option3_3Appear()
    AddXiuweiAppear()
    move1()
    move2()
    move3()
    move4()
    move5()
    move6()
    move7()
    move8()
    JingjieFormat()
    hint()
}, 50);
//自动保存(5秒一次)
window.setInterval(() => {
    save();
}, 100);
function upgrade1()
{
    if(userData.Lingqi.gte(Gongfa1cost[userData.Gongfa1Level]))userData.Lingqi = userData.Lingqi.sub(Gongfa1cost[userData.Gongfa1Level]),userData.Gongfa1Level++
    
}
function upgrade2()
{
    if(userData.Lingqi.gte(Gongfa2cost[userData.Gongfa2Level])&&userData.Jingjieid>=6)userData.Lingqi = userData.Lingqi.sub(Gongfa2cost[userData.Gongfa2Level]),userData.Gongfa2Level++
    
}
function upgrade3()
{
    if(userData.Lingqi.gte(Gongfa3cost[userData.Gongfa3Level])&&userData.Jingjieid>=11)userData.Lingqi = userData.Lingqi.sub(Gongfa3cost[userData.Gongfa3Level]),userData.Gongfa3Level++
    
}
function upgrade4()
{
    if(userData.Lingqi.gte(Gongfa4cost[userData.Gongfa4Level])&&userData.Jingjieid>=16)userData.Lingqi = userData.Lingqi.sub(Gongfa4cost[userData.Gongfa4Level]),userData.Gongfa4Level++
    
}
function upgrade5()
{
    if(userData.Lingqi.gte(Gongfa5cost[userData.Gongfa5Level])&&userData.Jingjieid>=21)userData.Lingqi = userData.Lingqi.sub(Gongfa5cost[userData.Gongfa5Level]),userData.Gongfa5Level++
    
}
function upgrade6()
{
    if(userData.Lingqi.gte(Gongfa6cost[userData.Gongfa6Level])&&userData.Jingjieid>=26)userData.Lingqi = userData.Lingqi.sub(Gongfa6cost[userData.Gongfa6Level]),userData.Gongfa6Level++
    
}
function upgrade7()
{
    if(userData.Lingqi.gte(Gongfa7cost[userData.Gongfa7Level])&&userData.Jingjieid>=31)userData.Lingqi = userData.Lingqi.sub(Gongfa7cost[userData.Gongfa7Level]),userData.Gongfa7Level++
    
}
function upgrade8()
{
    if(userData.Lingqi.gte(Gongfa8cost[userData.Gongfa8Level])&&userData.Jingjieid>=38)userData.Lingqi = userData.Lingqi.sub(Gongfa8cost[userData.Gongfa8Level]),userData.Gongfa8Level++
    
}
function upgrade9()
{
    if(userData.Lingqi.gte(Gongfa9cost[userData.Gongfa9Level])&&userData.Jingjieid>=46)userData.Lingqi = userData.Lingqi.sub(Gongfa9cost[userData.Gongfa9Level]),userData.Gongfa9Level++
    
}
function upgrade21()
{
    if(userData.Lingqi.gte(Julingzhen1cost[userData.Julingzhen1Level]))userData.Lingqi = userData.Lingqi.sub(Julingzhen1cost[userData.Julingzhen1Level]),userData.Julingzhen1Level++
    
}
function upgrade22()
{
    if(userData.Lingqi.gte(Julingzhen2cost[userData.Julingzhen2Level]))userData.Lingqi = userData.Lingqi.sub(Julingzhen2cost[userData.Julingzhen2Level]),userData.Julingzhen2Level++
    
}
function upgrade23()
{
    if(userData.Lingqi.gte(Julingzhen3cost[userData.Julingzhen3Level]))userData.Lingqi = userData.Lingqi.sub(Julingzhen3cost[userData.Julingzhen3Level]),userData.Julingzhen3Level++
    
}
function upgrade24()
{
    if(userData.Lingqi.gte(Julingzhen4cost[userData.Julingzhen4Level]))userData.Lingqi = userData.Lingqi.sub(Julingzhen4cost[userData.Julingzhen4Level]),userData.Julingzhen4Level++
    
}
function upgrade25()
{
    if(userData.Lingqi.gte(Julingzhen5cost[userData.Julingzhen5Level]))userData.Lingqi = userData.Lingqi.sub(Julingzhen5cost[userData.Julingzhen5Level]),userData.Julingzhen5Level++
    
}
function upgrade26()
{
    if(userData.Lingqi.gte(Julingzhen6cost[userData.Julingzhen6Level]))userData.Lingqi = userData.Lingqi.sub(Julingzhen6cost[userData.Julingzhen6Level]),userData.Julingzhen6Level++
    
}
function upgrade27()
{
    if(userData.Lingqi.gte(Julingzhen7cost[userData.Julingzhen7Level]))userData.Lingqi = userData.Lingqi.sub(Julingzhen7cost[userData.Julingzhen7Level]),userData.Julingzhen7Level++
    
}
function upgrade28()
{
    if(userData.Lingqi.gte(Julingzhen8cost[userData.Julingzhen8Level]))userData.Lingqi = userData.Lingqi.sub(Julingzhen8cost[userData.Julingzhen8Level]),userData.Julingzhen8Level++
    
}
function close0() {
    document.getElementById("option10").style.display = "none";
    document.getElementById("option20").style.display = "none";
    document.getElementById("option30").style.display = "none";
}
function close3() {
  document.getElementById("option3_1").style.display = "none";
  document.getElementById("option3_2").style.display = "none";
  document.getElementById("option3_3").style.display = "none";
}

function option1() {
    close0();
    document.getElementById("option10").style.display = "block";
}
function option2Appear() {
    if(userData.totalLingqi.gte(10000))document.getElementById("option2").style.display = "block";
    else document.getElementById("option2").style.display = "none";
}
function option3Appear() {
    if(userData.Lingshi.gte(150))document.getElementById("option3").style.display = "block";
    else document.getElementById("option3").style.display = "none";
}
function option4Appear() {
  if(userData.Jingjieid.gte(40))document.getElementById("option4").style.display = "block";
  else document.getElementById("option4").style.display = "none";
}
function AddXiuweiAppear() {
    if(userData.XiulianBalance.gt(1))document.getElementById("addXiuwei").style.display = "block";
    else document.getElementById("addXiuwei").style.display = "none";
}
function option3_2Appear() {
  if(userData.Jingjieid.gte(16))document.getElementById("option32").style.display = "block";
  else document.getElementById("option32").style.display = "none";
}
function option3_3Appear() {
  if(userData.Jingjieid.gte(31))document.getElementById("option33").style.display = "block";
  else document.getElementById("option33").style.display = "none";
}
function ZhujiAppear()
{
  if(userData.Jingjieid.gte(16))document.getElementById("Zhuji1").style.display = "block";
  else document.getElementById("Zhuji1").style.display = "none";
  if(userData.Jingjieid.gte(16))document.getElementById("Zhuji2").style.display = "block";
  else document.getElementById("Zhuji2").style.display = "none";
}
function KaiguangAppear()
{
  if(userData.Jingjieid.gte(31))document.getElementById("Kaiguang1").style.display = "block";
  else document.getElementById("Kaiguang1").style.display = "none";
  if(userData.Jingjieid.gte(31))document.getElementById("Kaiguang2").style.display = "block";
  else document.getElementById("Kaiguang2").style.display = "none";
}
function option2() {
    close0();
    document.getElementById("option20").style.display = "block";
}
function option3() {
    close0();
    document.getElementById("option30").style.display = "block";
}
function option3_1() {
  close3();
  document.getElementById("option3_1").style.display = "block";
}
function option3_2() {
  close3();
  document.getElementById("option3_2").style.display = "block";
}
function option3_3() {
  close3();
  document.getElementById("option3_3").style.display = "block";
}

function move1() {
    var elem = document.getElementById("myBar");   
        elem.style.width = format(format(userData.Julingzhen1FreezeLimit.sub(userData.Julingzhen1Freeze).div(userData.Julingzhen1FreezeLimit))*100) + '%'; 
        elem.innerHTML = format(format(userData.Julingzhen1FreezeLimit.sub(userData.Julingzhen1Freeze).div(userData.Julingzhen1FreezeLimit)) * 100)  + '%';
      }
      function move2() {
        var elem = document.getElementById("myBar2");   
            elem.style.width = format(format(userData.Julingzhen2FreezeLimit.sub(userData.Julingzhen2Freeze).div(userData.Julingzhen2FreezeLimit))*100) + '%'; 
            elem.innerHTML = format(format(userData.Julingzhen2FreezeLimit.sub(userData.Julingzhen2Freeze).div(userData.Julingzhen2FreezeLimit)) * 100)  + '%';
          }
          function move3() {
            var elem = document.getElementById("myBar3");   
                elem.style.width = format(format(userData.Julingzhen3FreezeLimit.sub(userData.Julingzhen3Freeze).div(userData.Julingzhen3FreezeLimit))*100) + '%'; 
                elem.innerHTML = format(format(userData.Julingzhen3FreezeLimit.sub(userData.Julingzhen3Freeze).div(userData.Julingzhen3FreezeLimit)) * 100)  + '%';
              }
              function move4() {
                var elem = document.getElementById("myBar4");   
                    elem.style.width = format(format(userData.Julingzhen4FreezeLimit.sub(userData.Julingzhen4Freeze).div(userData.Julingzhen4FreezeLimit))*100) + '%'; 
                    elem.innerHTML = format(format(userData.Julingzhen4FreezeLimit.sub(userData.Julingzhen4Freeze).div(userData.Julingzhen4FreezeLimit)) * 100)  + '%';
                    elem.style.backgroundColor = "#009600";
                  }
                  function move5() {
                    var elem = document.getElementById("myBar5");   
                        elem.style.width = format(format(userData.Julingzhen5FreezeLimit.sub(userData.Julingzhen5Freeze).div(userData.Julingzhen5FreezeLimit))*100) + '%'; 
                        elem.innerHTML = format(format(userData.Julingzhen5FreezeLimit.sub(userData.Julingzhen5Freeze).div(userData.Julingzhen5FreezeLimit)) * 100)  + '%';
                        elem.style.backgroundColor = "#009600";
                      }
                      function move6() {
                        var elem = document.getElementById("myBar6");   
                            elem.style.width = format(format(userData.Julingzhen6FreezeLimit.sub(userData.Julingzhen6Freeze).div(userData.Julingzhen6FreezeLimit))*100) + '%'; 
                            elem.innerHTML = format(format(userData.Julingzhen6FreezeLimit.sub(userData.Julingzhen6Freeze).div(userData.Julingzhen6FreezeLimit)) * 100)  + '%';
                            elem.style.backgroundColor = "#009600";
                          }
                          function move7() {
                            var elem = document.getElementById("myBar7");   
                                elem.style.width = format(format(userData.Julingzhen7FreezeLimit.sub(userData.Julingzhen7Freeze).div(userData.Julingzhen7FreezeLimit))*100) + '%'; 
                                elem.innerHTML = format(format(userData.Julingzhen7FreezeLimit.sub(userData.Julingzhen7Freeze).div(userData.Julingzhen7FreezeLimit)) * 100)  + '%';
                                elem.style.backgroundColor = "#0000FF";
                              }
                              function move8() {
                                var elem = document.getElementById("myBar8");   
                                    elem.style.width = format(format(userData.Julingzhen8FreezeLimit.sub(userData.Julingzhen8Freeze).div(userData.Julingzhen8FreezeLimit))*100) + '%'; 
                                    elem.innerHTML = format(format(userData.Julingzhen8FreezeLimit.sub(userData.Julingzhen8Freeze).div(userData.Julingzhen8FreezeLimit)) * 100)  + '%';
                                    elem.style.backgroundColor = "#0000FF";
                                  }
                      
  function JingjieFormat() {
    var elem = document.getElementById("XiuweiBar");   
        elem.style.width = format(userData.Xiuwei.div(JingjieLimit[userData.Jingjieid]))*70+30 + '%'; 
        if(parseInt((format(userData.Jingjieid)))<15.5) document.getElementById("XiuweiBar").style.backgroundColor = "#000000"
        else if(parseInt((format(userData.Jingjieid)))<30.5) document.getElementById("XiuweiBar").style.backgroundColor = "#009600"
        else if(parseInt((format(userData.Jingjieid)))<45.5) document.getElementById("XiuweiBar").style.backgroundColor = "#0000FF"
        elem.innerHTML = "当前境界："+userData.Jingjie+"(下阶所需"+format(userData.Xiuwei)+" / "+format(JingjieLimit[userData.Jingjieid])+"修为)";
        
    }

function hint()
{
    if(userData.totalLingqi.lt(10000)) document.getElementById("hint").innerHTML = "总共获得10000灵气以解锁聚灵阵。("+format(userData.totalLingqi.div(10000).mul(100).min(100))+"%)"
    else if(userData.Lingshi.lt(150))document.getElementById("hint").innerHTML = "到达150灵石以解锁练功房。("+format(userData.Lingshi.div(150).mul(100).min(100))+"%)"
    else if(userData.Jingjieid.lt(16))document.getElementById("hint").innerHTML = "境界到达筑基境可解锁更多练功房与功法选项。("+format(userData.Jingjieid.div(16).mul(100).min(100))+"%)"
    else if(userData.Jingjieid.lt(31))document.getElementById("hint").innerHTML = "境界到达开光境可解锁更多练功房与功法选项。("+format(userData.Jingjieid.div(31).mul(100).min(100))+"%)"
    else if(userData.Jingjieid.lt(40))document.getElementById("hint").innerHTML = "境界到达开光10阶可解锁转生。("+format(userData.Jingjieid.div(40).mul(100).min(100))+"%)"
    else document.getElementById("hint").innerHTML = "当前残局：100亿灵气&&境界达到开光10阶（第一次转生要求）"
}
function Xiulian11()
{
    if(userData.Lingshi.gte(150))userData.XiulianBalance = new Decimal(15),userData.XiuweiPerSec = new Decimal(100),userData.Lingshi = userData.Lingshi.sub(150)
}
function Xiulian12()
{
    if(userData.Lingshi.gte(245))userData.XiulianBalance = new Decimal(15),userData.XiuweiPerSec = new Decimal(125),userData.Lingshi = userData.Lingshi.sub(245)
}
function Xiulian13()
{
    if(userData.Lingshi.gte(900))userData.XiulianBalance = new Decimal(15),userData.XiuweiPerSec = new Decimal(225),userData.Lingshi = userData.Lingshi.sub(900)
}
function Xiulian14()
{
    if(userData.Lingshi.gte(5500))userData.XiulianBalance = new Decimal(15),userData.XiuweiPerSec = new Decimal(600),userData.Lingshi = userData.Lingshi.sub(5500)
    
    
}
function Xiulian21()
{
    if(userData.Lingshi.gte(15000))userData.XiulianBalance = new Decimal(15),userData.XiuweiPerSec = new Decimal(2000),userData.Lingshi = userData.Lingshi.sub(15000)
}
function Xiulian22()
{
    if(userData.Lingshi.gte(33000))userData.XiulianBalance = new Decimal(15),userData.XiuweiPerSec = new Decimal(3200),userData.Lingshi = userData.Lingshi.sub(33000)
}
function Xiulian23()
{
    if(userData.Lingshi.gte(90000))userData.XiulianBalance = new Decimal(15),userData.XiuweiPerSec = new Decimal(4900),userData.Lingshi = userData.Lingshi.sub(90000)
}
function Xiulian24()
{
    if(userData.Lingshi.gte(650000))userData.XiulianBalance = new Decimal(15),userData.XiuweiPerSec = new Decimal(13500),userData.Lingshi = userData.Lingshi.sub(650000)
}
function Xiulian31()
{
    if(userData.Lingshi.gte(420000))userData.XiulianBalance = new Decimal(15),userData.XiuweiPerSec = new Decimal(27000),userData.Lingshi = userData.Lingshi.sub(420000)
}
function Xiulian32()
{
    if(userData.Lingshi.gte(880000))userData.XiulianBalance = new Decimal(15),userData.XiuweiPerSec = new Decimal(44000),userData.Lingshi = userData.Lingshi.sub(880000)
}
function Xiulian33()
{
    if(userData.Lingshi.gte(2880000))userData.XiulianBalance = new Decimal(15),userData.XiuweiPerSec = new Decimal(120000),userData.Lingshi = userData.Lingshi.sub(2880000)
}
function Xiulian34()
{
    if(userData.Lingshi.gte(10000000))userData.XiulianBalance = new Decimal(15),userData.XiuweiPerSec = new Decimal(300000),userData.Lingshi = userData.Lingshi.sub(10000000)
}

function getAchievement() {
  document.getElementById("getAchievement").style.cssText = "opacity: 0.1";
    document.getElementById("getAchievement").style.display = "block";
    setTimeout("document.getElementById('getAchievement').style.cssText = 'opacity: 0.2';",50)
    setTimeout("document.getElementById('getAchievement').style.cssText = 'opacity: 0.3';",100)
    setTimeout("document.getElementById('getAchievement').style.cssText = 'opacity: 0.4';",150)
    setTimeout("document.getElementById('getAchievement').style.cssText = 'opacity: 0.5';",200)
    setTimeout("document.getElementById('getAchievement').style.cssText = 'opacity: 0.6';",250)
    setTimeout("document.getElementById('getAchievement').style.cssText = 'opacity: 0.7';",300)
    setTimeout("document.getElementById('getAchievement').style.cssText = 'opacity: 0.8';",350)
    setTimeout("document.getElementById('getAchievement').style.cssText = 'opacity: 0.9';",400)
    setTimeout("document.getElementById('getAchievement').style.cssText = 'opacity: 1.0';",450)
    setTimeout("document.getElementById('getAchievement').style.cssText = 'opacity: 0.9';",2200)
    setTimeout("document.getElementById('getAchievement').style.cssText = 'opacity: 0.8';",2250)
    setTimeout("document.getElementById('getAchievement').style.cssText = 'opacity: 0.7';",2300)
    setTimeout("document.getElementById('getAchievement').style.cssText = 'opacity: 0.6';",2350)
    setTimeout("document.getElementById('getAchievement').style.cssText = 'opacity: 0.5';",2400)
    setTimeout("document.getElementById('getAchievement').style.cssText = 'opacity: 0.4';",2450)
    setTimeout("document.getElementById('getAchievement').style.cssText = 'opacity: 0.3';",2500)
    setTimeout("document.getElementById('getAchievement').style.cssText = 'opacity: 0.2';",2550)
    setTimeout("document.getElementById('getAchievement').style.cssText = 'opacity: 0.1';",2600)
    setTimeout("document.getElementById('getAchievement').style.display = 'none';",2650)
}