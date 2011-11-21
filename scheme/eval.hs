module Eval (eval)
where

import Types
import Dispatch
import Control.Monad.Error
import Data.Map (empty)

eval :: LispExpr -> Env -> LispMonad LispExpr
eval expr env = case expr of
  Primitive _   -> return expr
  Closure _ _ _ -> return expr
  Atom sym      -> lookupEnv sym env
  List l        -> listEval l env
  _             -> malformed expr

listEval :: [LispExpr] -> Env -> LispMonad LispExpr
listEval [] _ = return $ List []

listEval l@(Atom h:t) env =
  if reserved h
  then dispatch h eval env (List l)
  else do head <- eval (Atom h) env
          apply head t env

listEval (h:t) env =
  do head <- eval h env
     apply head t env

apply :: LispExpr -> [LispExpr] -> Env -> LispMonad LispExpr
apply expr [] _ = return expr

apply (Closure expr farg env') (arg:args) env =
  do argv <- eval arg env
     closure <- eval expr (bind farg argv env')
     apply closure args env

apply expr args _ = throwError $ NonFunctional expr args