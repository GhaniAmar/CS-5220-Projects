module Dispatch
where

import Keywords
import Operators
import Data.Map (fromList, member, (!))

keywordTable =
  fromList  [("apply", applyKeyword),
             ("car", carKeyword),
             ("cdr", cdrKeyword),
             ("cond", condKeyword),
             ("cons", consKeyword),
             ("display", displayKeyword),
             ("eq?", eqKeyword),
             ("equal?", equalKeyword),
             ("eval", evalKeyword),
             ("if", ifKeyword),
             ("lambda", lambdaKeyword),
             ("let", letKeyword),
             ("letrec", letrecKeyword),
             ("list", listKeyword),
             ("list?", listqKeyword),
             ("null?", nullKeyword),
             ("parallel-let", pletKeyword),
             ("parallel-map", pMapKeyword),
             ("quote", quoteKeyword)]

operatorTable =
  fromList [("*", timesOperator),
            ("+", plusOperator),
            ("-", minusOperator),
            ("/", divideOperator),
            ("<", lessOperator),
            ("<=", leqOperator),
            ("=", eqOperator),
            (">", gtrOperator),
            (">=", geqOperator),
            ("and", andOperator),
            ("not", notOperator),
            ("or", orOperator)]

reserved s = member s keywordTable || member s operatorTable
dispatch s | member s keywordTable = keywordTable ! s
           | otherwise             = operatorTable ! s