module Operators
where

import Types
import Control.Monad.Error

data TypeResult = Int | Flt | Fail

checkType :: [LispExpr] -> TypeResult
checkType = foldl helper Int
  where helper Flt (Primitive (Integer _)) = Flt
        helper Flt (Primitive (Double _))  = Flt
        helper Int (Primitive (Integer _)) = Int
        helper Int (Primitive (Double _))  = Flt
        helper _ _                         = Fail

toFlt :: LispExpr -> LispMonad Double
toFlt (Primitive (Integer n)) = return . fromInteger $ n
toFlt (Primitive (Double n))  = return n
toFlt e                       = throwError . TypeError $ e

toInt :: LispExpr -> LispMonad Integer
toInt (Primitive (Integer e)) = return e
toInt (Primitive (Double e))  = return . floor $ e
toInt e                       = throwError . TypeError $ e

numericFunction op1 op2 check eval env l@(List (Atom _ : t)) =
  do args <- mapM (`eval` env) t
     verify args
     case checkType args of
       Int ->
         do ints <- mapM toInt args
            return . Primitive . Integer $ foldl1 op1 ints
       Flt ->
         do floats <- mapM toFlt args
            return . Primitive . Double $ foldl1 op2 floats
       Fail -> throwError . TypeError $ l
  where verify ls = case ls of
          x:y:_ -> mapM (check l) (tail ls)
          _     -> return []

numericFunction _ _ _ _ _ l = malformed l

boolFunction op eval env (List (Atom _ : t)) =
  do args <- mapM (`eval` env) t
     let argTypes = checkType args
     coerced <- mapM (coerce argTypes) args
     return . Primitive . Bool $ case coerced of
       x:xs -> fst (foldl folder (True, x) xs)
       []   -> True
  where folder (b, lst) e = (op b lst e, e)
        coerce Flt e = do x <- toFlt e
                          return . Primitive . Double $ x
        coerce _ x   = return x
boolFunction _ _ _ l = malformed l

wrap op b lst e = b && op lst e
vacuous _ _ = return True
byZero e (Primitive (Integer 0)) = throwError . NumericError $ e
byZero e (Primitive (Double 0))  = throwError . NumericError $ e
byZero e _ = return True

timesOperator  = numericFunction (*) (*) vacuous
plusOperator   = numericFunction (+) (+) vacuous
minusOperator  = numericFunction (-) (-) vacuous
divideOperator = numericFunction div (/) byZero
lessOperator   = boolFunction (wrap (<))
leqOperator    = boolFunction (wrap (<=))
eqOperator     = boolFunction (wrap (==))
gtrOperator    = boolFunction (wrap (>))
geqOperator    = boolFunction (wrap (>=))

andOperator eval env (List (Atom "and" : args)) =
  do args' <- mapM (\x -> eval x env >>= toBool) args
     return . Primitive . Bool . and $ args'
andOperator _ _ l = malformed l

orOperator eval env (List (Atom "or" : args)) =
  do args' <- mapM (\x -> eval x env >>= toBool) args
     return . Primitive . Bool . or $ args'
orOperator _ _ l = malformed l

notOperator eval env (List [Atom "not", e]) =
  do result <- toBool =<< eval e env
     return . Primitive . Bool . not $ result
notOperator _ _ l = malformed l